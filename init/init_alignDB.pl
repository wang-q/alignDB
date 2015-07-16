#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use Text::CSV_XS;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB;
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

# Database init values
my $server   = $Config->{database}{server};
my $port     = $Config->{database}{port};
my $username = $Config->{database}{username};
my $password = $Config->{database}{password};
my $db       = $Config->{database}{db};

my $init_sql   = "$FindBin::Bin/../init.sql";
my $init_taxon = "$FindBin::Bin/../data/taxon.csv";
my $init_chr   = "$FindBin::Bin/../data/chr_length.csv";

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'             => \$help,
    'man'                => \$man,
    's|server=s'         => \$server,
    'P|port=i'           => \$port,
    'd|db=s'             => \$db,
    'u|username=s'       => \$username,
    'p|password=s'       => \$password,
    'sql|init_sql=s'     => \$init_sql,
    'taxon|init_taxon=s' => \$init_taxon,
    'chr|init_chr=s'     => \$init_chr,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# call mysql
#----------------------------------------------------------#
{
    $stopwatch->block_message("Create DB skeleton");

    my $drh = DBI->install_driver("mysql");    # Driver handle object
    $drh->func( 'dropdb',   $db, $server, $username, $password, 'admin' );
    $drh->func( 'createdb', $db, $server, $username, $password, 'admin' );

    my $dbh = DBI->connect( "dbi:mysql:$db:$server", $username, $password );
    open my $infh, '<', $init_sql;
    my $content = do { local $/; <$infh> };
    close $infh;
    my @statements = grep {/\w/} split /;/, $content;
    for (@statements) {
        $dbh->do($_) or die $dbh->errstr;
    }
}

#----------------------------------------------------------#
# init object
#----------------------------------------------------------#
$stopwatch->start_message("Init $db...");

my $obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);
$obj->index_isw_indel_id;
my $dbh = $obj->dbh;

#----------------------------------------------------------#
# taxon
#----------------------------------------------------------#
if ( -f $init_taxon ) {
    print "Use $init_taxon to Init table taxon\n";

    my $insert_sth = $dbh->prepare(
        'INSERT INTO taxon (
            taxon_id, genus, species, sub_species, common_name
        )
        VALUES (
            ?, ?, ?, ?, ? 
        )'
    );

    my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n" } );
    open my $csv_fh, "<", $init_taxon;
    $csv->getline($csv_fh);    # bypass title line

    while ( my $row = $csv->getline($csv_fh) ) {
        $insert_sth->execute( $row->[0], $row->[1], $row->[2], $row->[3],
            $row->[4], );
    }
    close $csv_fh;
    $insert_sth->finish;

    # Now retrieve data from the table.
    my $query_sth = $dbh->prepare(
        q{
        SELECT *
        FROM taxon
        WHERE common_name = "Human"
        }
    );
    $query_sth->execute;
    while ( my $ref = $query_sth->fetchrow_hashref ) {
        printf "Found a row: taxon_id = %s, genus = %s,"
            . " species = %s, common_name = %s\n",
            $ref->{taxon_id},
            $ref->{genus}, $ref->{species}, $ref->{common_name};
    }
    $query_sth->finish;

    print "\n";
}

#----------------------------------------------------------#
# chromosome
#----------------------------------------------------------#
if ( -f $init_chr ) {
    print "Use $init_chr to Init table chromosome\n";

    my $insert_sth = $dbh->prepare(
        'INSERT INTO chromosome (
            chr_id, taxon_id, chr_name, chr_length
        )
        VALUES (
            NULL, ?, ?, ?
        )'
    );

    my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n" } );
    open my $csv_fh, "<", $init_chr;
    $csv->getline($csv_fh);    # bypass title line

    while ( my $row = $csv->getline($csv_fh) ) {
        $insert_sth->execute( $row->[0], $row->[1], $row->[2], );
    }
    close $csv_fh;
    $insert_sth->finish;

    # Now retrieve data from the table
    my $query_sth = $dbh->prepare(
        q{
        SELECT *
        FROM chromosome
        }
    );
    $query_sth->execute;
    while ( my $ref = $query_sth->fetchrow_hashref ) {
        my $taxon_id = $ref->{taxon_id};
        if ( defined $taxon_id and $taxon_id == 3702 ) {
            printf "Found a row: taxon_id = %s,"
                . " chr_name = %s, chr_length = %s\n",
                $ref->{taxon_id}, $ref->{chr_name}, $ref->{chr_length};
        }
    }
    $query_sth->finish;
}

$stopwatch->end_message;

# store program running meta info to database
# this AlignDB object is just for storing meta info
END {
    AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    )->add_meta_stopwatch($stopwatch);
}

__END__


=head1 NAME

init_alignDB.pl - Initiate alignDB

=head1 SYNOPSIS

    perl init_alignDB.pl [options]
      Options:
        --help          brief help message
        --man           full documentation
        --server        MySQL server IP/Domain name
        --port          MySQL server port
        --db            database name
        --username      username
        --password      password
        --init_sql      init sql filename

    perl init_alignDB.pl -d S288cvsYJM789

=cut
