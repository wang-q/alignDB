#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use Path::Tiny;
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
    'help|?' => sub { Getopt::Long::HelpMessage(0) },
    'server|s=s'   => \( my $server   = $Config->{database}{server} ),
    'port=i'       => \( my $port     = $Config->{database}{port} ),
    'db|d=s'       => \( my $db       = $Config->{database}{db} ),
    'username|u=s' => \( my $username = $Config->{database}{username} ),
    'password|p=s' => \( my $password = $Config->{database}{password} ),
    'sql=s'        => \( my $init_sql = "$FindBin::RealBin/../init.sql" ),
    'chr=s'        => \( my $init_chr = $Config->{generate}{file_chr_length} ),
) or Getopt::Long::HelpMessage(1);

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
    my DBI $dbh = DBI->connect( "dbi:mysql:$db:$server", $username, $password );
    my $in_fh = path($init_sql)->openr;
    my $content = do { local $/; <$in_fh> };
    close $in_fh;
    my @statements = grep {/\w/} split /;/, $content;
    for (@statements) {
        $dbh->do($_) or die $dbh->errstr;
    }
}

#----------------------------------------------------------#
# chromosome
#----------------------------------------------------------#
{
    $stopwatch->block_message("Use [$init_chr] to Init table chromosome");

    my $obj = AlignDB->new(
        mysql  => "$db:$server",
        user   => $username,
        passwd => $password,
    );
    my DBI $dbh = $obj->dbh;

    my DBI $insert_sth = $dbh->prepare(
        'INSERT INTO chromosome (
            chr_id, common_name, taxon_id, chr_name, chr_length
        )
        VALUES (
            NULL, ?, ?, ?, ?
        )'
    );

    my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n" } );
    my $csv_fh = path($init_chr)->openr;
    $csv->getline($csv_fh);    # bypass title line

    while ( my $row = $csv->getline($csv_fh) ) {
        $insert_sth->execute( $row->[0], $row->[1], $row->[2], $row->[3], );
    }
    close $csv_fh;
    $insert_sth->finish;
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
