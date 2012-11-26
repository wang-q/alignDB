#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use DBI;

use AlignDB::Stopwatch;

use FindBin;

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

# Database init values
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db       = $Config->{bac}{db};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'          => \$help,
    'man'             => \$man,
    'server=s'        => \$server,
    'port=i'          => \$port,
    'db=s'            => \$db,
    'username=s'      => \$username,
    'password=s'      => \$password,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Blastz whole species");

my $dbh = DBI->connect( "dbi:mysql:$db:$server", $username, $password );

my @taxon_ids;
{    # taxon.csv

    # select all strains in this species
    my $sth = $dbh->prepare(
        qq{
        SELECT  st.taxonomy_id,
                st.genus,
                st.species,
                st.organism_name,
                st.released_date
        FROM strain st
        WHERE st.seq_ok = 1
        ORDER BY st.group,
                st.species_id,
                st.released_date
        }
    );
    $sth->execute;

    # "taxon_id,genus,species,sub_species,common_name,classification\n";
    open my $fh, '>', "taxon.csv";
    while ( my @row = $sth->fetchrow_array ) {
        my $taxon_id = $row[0];
        my $genus = $row[1];
        my $species = $row[2];
        my $sub_species = $row[3];
        
        $sub_species =~ s/^$species//;
        $sub_species =~ s/^\s+//;
        
        $species =~ s/^$genus//;
        $species =~ s/^\s+//;
        
        print {$fh} "$taxon_id,$genus,$species,$sub_species\n";
        
        push @taxon_ids, $taxon_id;
    }
    close $fh;
}

{    # chr_length.csv

    # "taxon_id,chr,length,name,assembly\n";
    open my $fh, '>', "chr_length.csv";
    
    for my $taxon_id ( @taxon_ids ) {

        my $sth = $dbh->prepare(
            qq{
            SELECT  s.accession,
                    s.organism_name,
                    s.replicon,
                    s.length
            FROM seq s
            WHERE s.taxonomy_id = ?
            }
        );
        $sth->execute($taxon_id);

        while ( my ( $acc, $name, $rep, $length ) = $sth->fetchrow_array ) {

            # accession as chr name
            # taxon id as strain name
            # replicon as assembly name
            $name =~ s/\W/_/g;
            $rep =~ s/^chromosome/chr/;
            $rep =~ s/^plasmid p/p/;
            $rep =~ s/^plasmid/p/;
            $rep =~ s/\s+//g;
            print {$fh}
                "$taxon_id,$acc,$length,${name}_${taxon_id},$rep\n";
        }
    }
    close $fh;
}

$stopwatch->end_message;
exit;

__END__

perl bac_pre_aligndb.pl 
