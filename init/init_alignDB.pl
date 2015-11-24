#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use Text::CSV_XS;

use lib "$FindBin::Bin/../lib";
use AlignDB;
use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# record ARGV and Config
my $stopwatch = AlignDB::Stopwatch->new(
    program_name => $0,
    program_argv => [@ARGV],
    program_conf => $Config,
);

=head1 NAME

init_alignDB.pl - Initiate alignDB

=head1 SYNOPSIS

    perl init_alignDB.pl [options]
      Options:
        --help      -?          brief help message
        --server    -s  STR     MySQL server IP/Domain name
        --port      -P  INT     MySQL server port
        --db        -d  STR     database name
        --username  -u  STR     username
        --password  -p  STR     password
        --sql           STR     init sql filename
        --taxon         STR     init taxon filename
        --chr           STR     init chr_length filename

    perl init_alignDB.pl -d S288cvsYJM789

=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'server|s=s'   => \( my $server     = $Config->{database}{server} ),
    'port|P=i'     => \( my $port       = $Config->{database}{port} ),
    'db|d=s'       => \( my $db         = $Config->{database}{db} ),
    'username|u=s' => \( my $username   = $Config->{database}{username} ),
    'password|p=s' => \( my $password   = $Config->{database}{password} ),
    'sql=s'        => \( my $init_sql   = "$FindBin::Bin/../init.sql" ),
    'taxon=s'      => \( my $init_taxon = "$FindBin::Bin/../data/taxon.csv" ),
    'chr=s'        => \( my $init_chr   = "$FindBin::Bin/../data/chr_length.csv" ),
) or HelpMessage(1);

#----------------------------------------------------------#
# call mysql
#----------------------------------------------------------#
{
    $stopwatch->block_message("Create DB skeleton");

    my $drh = DBI->install_driver("mysql");    # Driver handle object
    print "# dropdb\n";
    $drh->func( 'dropdb', $db, $server, $username, $password, 'admin' );
    print "# createdb\n";
    $drh->func( 'createdb', $db, $server, $username, $password, 'admin' );

    print "# init\n";
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
        $insert_sth->execute( $row->[0], $row->[1], $row->[2], $row->[3], $row->[4], );
    }
    close $csv_fh;
    $insert_sth->finish;

    # Now retrieve data from the table.
    my $query_sth = $dbh->prepare(
        q{
        SELECT *
        FROM taxon
        WHERE common_name in ( "Human", "human", "S288C", "S288c")
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
        if ( defined $taxon_id and ( $taxon_id == 3702 or $taxon_id == 9606 ) ) {
            printf "Found a row: taxon_id = %s,"
                . " chr_name = %s, chr_length = %s\n",
                $ref->{taxon_id}, $ref->{chr_name}, $ref->{chr_length};
            last;
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
