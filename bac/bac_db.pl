#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use DBI;
use Text::CSV_XS;

use FindBin;

use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

# running options
my $strain_file  = $Config->{bac}{strain_file};
my $seq_file     = $Config->{bac}{seq_file};
my $species_file = $Config->{bac}{species_file};

# Database init values
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db       = $Config->{bac}{db};

my $init_sql = "$FindBin::Bin/../bac.sql";

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'     => \$help,
    'man'        => \$man,
    'server=s'   => \$server,
    'port=i'     => \$port,
    'db=s'       => \$db,
    'username=s' => \$username,
    'password=s' => \$password,
    'init_sql=s' => \$init_sql,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Init bac DB...");

#----------------------------#
# call mysql
#----------------------------#
{
    $stopwatch->block_message("Create DB skeleton");

    my $drh = DBI->install_driver("mysql");    # Driver handle object
    $drh->func( 'dropdb',   $db, $server, $username, $password, 'admin' );
    $drh->func( 'createdb', $db, $server, $username, $password, 'admin' );
    $drh->func( 'reload',   $db, $server, $username, $password, 'admin' );

    my $dbh = DBI->connect( "dbi:mysql:$db:$server", $username, $password );
    open my $infh, '<', $init_sql;
    my $content = do { local $/; <$infh> };
    close $infh;
    my @statements = grep {/\w/} split /;/, $content;
    for (@statements) {
        $dbh->do($_) or die $dbh->errstr;
    }
}

#----------------------------#
# Filling table strain
#----------------------------#
my $dbh = DBI->connect( "dbi:mysql:$db:$server", $username, $password );

{
    $stopwatch->block_message("Loading $strain_file");

    my $load_sth = $dbh->prepare(
        qq{
        LOAD DATA LOCAL INFILE '$strain_file'
        INTO TABLE strain
        FIELDS TERMINATED BY ',' ENCLOSED BY '"' LINES TERMINATED BY '\n'
        IGNORE 1 LINES
        }
    );
    $load_sth->execute;
    $load_sth->finish;
}

{
    $stopwatch->block_message("Update species and genus memberships");

    # find species contains multiply strains
    my %species_member_of;
    my $species_sth = $dbh->prepare(
        qq{
        SELECT species, count(taxonomy_id)
        FROM   strain
        GROUP BY species
        }
    );
    $species_sth->execute;
    while ( my @row = $species_sth->fetchrow_array ) {
        my ( $species, $count ) = @row;
        $species_member_of{$species} = $count;
    }
    $species_sth->finish;

    # find genus contains multiply species
    my %genus_member_of;
    my $genus_sth = $dbh->prepare(
        qq{
        SELECT genus, count(distinct species), count(taxonomy_id)
        FROM   strain
        GROUP BY genus
        }
    );
    $genus_sth->execute;
    while ( my @row = $genus_sth->fetchrow_array ) {
        my ( $genus, $species_count, $strain_count ) = @row;
        $genus_member_of{$genus} = [ $species_count, $strain_count ];
    }
    $genus_sth->finish;

    my $update_sth = $dbh->prepare(
        qq{
        UPDATE  strain
        SET     species_member = ?,
                genus_species_member = ?,
                genus_strain_member = ?
        WHERE   taxonomy_id = ?
        }
    );
    my $id_sth = $dbh->prepare(
        qq{
        SELECT taxonomy_id, species, genus
        FROM   strain
        }
    );
    $id_sth->execute;
    while ( my @row = $id_sth->fetchrow_array ) {
        my ( $taxonomy_id, $species, $genus ) = @row;
        $update_sth->execute(
            $species_member_of{$species},  $genus_member_of{$genus}->[0],
            $genus_member_of{$genus}->[1], $taxonomy_id
        );
    }
    $id_sth->finish;
    $update_sth->finish;
}

#----------------------------#
# Filling table seq
#----------------------------#
{
    $stopwatch->block_message("Loading $seq_file");

    my $load_sth = $dbh->prepare(
        qq{
        LOAD DATA LOCAL INFILE '$seq_file'
        INTO TABLE seq
        FIELDS TERMINATED BY ',' ENCLOSED BY '"' LINES TERMINATED BY '\n'
        IGNORE 1 LINES
        }
    );
    $load_sth->execute;
    $load_sth->finish;
}

#----------------------------#
# Output a summary csv file
#----------------------------#
{
    $stopwatch->block_message("Writing $species_file");

    # prepare output csv file
    my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n" } );
    open my $csv_fh, ">", $species_file or die "$species_file: $!";
    my @headers
        = qw{species chr_number avg_genome_size avg_gc species_member genus_species_member genus_strain_member};
    $csv->print( $csv_fh, \@headers );

    # print species' member, chr_number and genus_member
    my $species_count_sth = $dbh->prepare(
        qq{
        SELECT  species,
                number_of_chromosomes,
                avg(genome_size),
                avg(gc_content),
                species_member,
                genus_species_member,
                genus_strain_member
        FROM    strain
        WHERE   number_of_chromosomes > 0
        GROUP BY species
        }
    );
    $species_count_sth->execute;
    while ( my @row = $species_count_sth->fetchrow_array ) {
        $csv->print( $csv_fh, \@row );
    }
    $species_count_sth->finish;
    close $csv_fh;
}

#----------------------------#
# check if strain matches seq
#----------------------------#
{
    $stopwatch->block_message("Checking refseq_accessions");

    my @nok_ids;    # wrong taxonomy_ids

    my $strain_sth = $dbh->prepare(
        qq{
        SELECT s.taxonomy_id id, s.organism_name name, s.refseq_accessions acc
        FROM strain s
        }
    );
    $strain_sth->execute;
    my $strain_ref = $strain_sth->fetchall_hashref('id');
    $strain_sth->finish;

    my $seq_sth = $dbh->prepare(
        qq{
        SELECT accession acc
        from seq
        where taxonomy_id = ?
        order by acc
        }
    );
    for my $id ( keys %{$strain_ref} ) {
        my $name = $strain_ref->{$id}{name};
        my $strain_str = join ",", sort split /,/, $strain_ref->{$id}{acc};

        $seq_sth->execute($id);
        my $ary_ref  = $seq_sth->fetchall_arrayref;
        my @seq_accs = map { ( @{$_} ) } @{$ary_ref};
        my $seq_str  = join ",", sort @seq_accs;

        if ( $strain_str ne $seq_str ) {
            warn "For taxonomy_id $id,  strain doesn't match seq\n";
            warn Dump {
                name__ => $name,
                seq___ => $seq_str,
                strain => $strain_str,
            };
            push @nok_ids, $id;
        }
    }
    $seq_sth->finish;

    $dbh->do(
        qq{
        ALTER TABLE strain
        Add COLUMN seq_ok int default 1
        }
    );
    my $update_sth = $dbh->prepare(
        qq{
        update strain
        set seq_ok = 0
        where taxonomy_id = ?
        }
    );
    for (@nok_ids) {
        $update_sth->execute($_);
    }
    $update_sth->finish;
}

#----------------------------#
# Finish
#----------------------------#
$stopwatch->end_message;
exit;

__END__

perl bac_db.pl 
