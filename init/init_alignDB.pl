#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use Text::CSV_XS;

use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/../lib";
use AlignDB;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini");

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
        --chr           STR     init chr_length filename

    perl init_alignDB.pl -d S288cvsYJM789

=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'server|s=s'   => \( my $server   = $Config->{database}{server} ),
    'port|P=i'     => \( my $port     = $Config->{database}{port} ),
    'db|d=s'       => \( my $db       = $Config->{database}{db} ),
    'username|u=s' => \( my $username = $Config->{database}{username} ),
    'password|p=s' => \( my $password = $Config->{database}{password} ),
    'sql=s'        => \( my $init_sql = "$FindBin::RealBin/../init.sql" ),
    'chr=s'        => \( my $init_chr = "$FindBin::RealBin/../data/chr_length.csv" ),
) or HelpMessage(1);

#----------------------------------------------------------#
# call mysql
#----------------------------------------------------------#
$stopwatch->start_message("Init $db...");

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
my $obj = AlignDB->new(
    mysql  => "$db:$server",
    user   => $username,
    passwd => $password,
);
my $dbh = $obj->dbh;

#----------------------------------------------------------#
# chromosome
#----------------------------------------------------------#
if ( -f $init_chr ) {
    $stopwatch->block_message("Use [$init_chr] to Init table chromosome");

    my $insert_sth = $dbh->prepare(
        'INSERT INTO chromosome (
            chr_id, common_name, taxon_id, chr_name, chr_length
        )
        VALUES (
            NULL, ?, ?, ?, ?
        )'
    );

    my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n" } );
    open my $csv_fh, "<", $init_chr;
    $csv->getline($csv_fh);    # bypass title line

    while ( my $row = $csv->getline($csv_fh) ) {
        $insert_sth->execute( $row->[0], $row->[1], $row->[2], $row->[3], );
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
        my $common_name = $ref->{common_name};
        if ( defined $common_name and ( $common_name eq "Human" ) ) {
            printf "Found a row: common_name = %s, chr_name = %s, chr_length = %s\n",
                $ref->{common_name}, $ref->{chr_name}, $ref->{chr_length};
            last;
        }
    }
    $query_sth->finish;
}
else {
    die "[$init_chr] doesn't exist\n";
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
