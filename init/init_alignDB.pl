#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long::Descriptive;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use Path::Tiny;
use Text::CSV_XS;

use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/../lib";
use AlignDB::Common;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $conf = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini");

# record command line
my $stopwatch = AlignDB::Stopwatch->new->record;

my $description = <<'EOF';
Initiate alignDB

    perl init/init_alignDB.pl -d S288CvsRM11_1a

Usage: perl %c [options]
EOF

(
    #@type Getopt::Long::Descriptive::Opts
    my $opt,

    #@type Getopt::Long::Descriptive::Usage
    my $usage,
    )
    = Getopt::Long::Descriptive::describe_options(
    $description,
    [ 'help|h', 'display this message' ],
    [],
    ['Database init values'],
    [ 'server|s=s',   'MySQL IP/Domain', { default => $conf->{database}{server} }, ],
    [ 'port=i',       'MySQL port',      { default => $conf->{database}{port} }, ],
    [ 'username|u=s', 'username',        { default => $conf->{database}{username} }, ],
    [ 'password|p=s', 'password',        { default => $conf->{database}{password} }, ],
    [ 'db|d=s',       'database name',   { default => $conf->{database}{db} }, ],
    [],
    [ 'sql=s', 'init sql filename',        { default => "$FindBin::RealBin/../init.sql" }, ],
    [ 'chr=s', 'init chr_length filename', { default => $conf->{generate}{file_chr_length} }, ],
    { show_defaults => 1, }
    );

$usage->die if $opt->{help};

# record config
$stopwatch->record_conf($opt);

# DBI Data Source Name
my $dsn = sprintf "dbi:mysql:database=%s;host=%s;port=%s", $opt->{db}, $opt->{server}, $opt->{port};

#----------------------------------------------------------#
# call mysql
#----------------------------------------------------------#
$stopwatch->start_message("Init [$opt->{db}]...");

{
    $stopwatch->block_message("Create DB skeleton");

    $ENV{MYSQL_PWD} = $opt->{password};
    my $cmd = "mysql -h$opt->{server} -P$opt->{port} -u$opt->{username}";

    my $drop = " -e \"DROP DATABASE IF EXISTS $opt->{db};\"";
    print "* dropdb\n" . "$cmd $drop\n";
    system("$cmd $drop");

    my $create = " -e \"CREATE DATABASE $opt->{db};\"";
    print "* createdb\n" . "$cmd $create\n";
    system("$cmd $create");

    #@type DBI
    my $dbh = DBI->connect( $dsn, $opt->{username}, $opt->{password}, )
        or Carp::confess $DBI::errstr;

    print "* init\n";
    my $content = path( $opt->{sql} )->slurp;
    my @statements = grep {/\w/} split /;/, $content;
    for (@statements) {
        $dbh->do($_) or die $dbh->errstr;
    }
}

#----------------------------------------------------------#
# chromosome
#----------------------------------------------------------#
{
    $stopwatch->block_message("Use [$opt->{chr}] to Init table chromosome");

    my $obj = AlignDB::Common->new(
        dsn    => $dsn,
        user   => $opt->{username},
        passwd => $opt->{password},
    );

    #@type DBI
    my $sth = $obj->dbh->prepare(
        'INSERT INTO chromosome (
            chr_id, common_name, taxon_id, chr_name, chr_length
        )
        VALUES (
            NULL, ?, ?, ?, ?
        )'
    );

    my $csv = Text::CSV_XS->new( { binary => 1, eol => "\n" } );
    my $csv_fh = path( $opt->{chr} )->openr;
    $csv->getline($csv_fh);    # bypass title line

    while ( my $row = $csv->getline($csv_fh) ) {
        $sth->execute( $row->[0], $row->[1], $row->[2], $row->[3], );
    }
    close $csv_fh;
    $sth->finish;
}

$stopwatch->end_message;

# store program's meta info to database
AlignDB::Common->new(
    dsn    => $dsn,
    user   => $opt->{username},
    passwd => $opt->{password},
)->add_meta_stopwatch($stopwatch);

__END__
