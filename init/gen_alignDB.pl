#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long::Descriptive;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use File::Find::Rule;
use MCE;
use Path::Tiny;

use AlignDB::Stopwatch;

use lib "$FindBin::RealBin/../lib";
use AlignDB;
use AlignDB::Outgroup;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $conf = Config::Tiny->read("$FindBin::RealBin/../alignDB.ini");

# record command line
my $stopwatch = AlignDB::Stopwatch->new->record;

my $description = <<'EOF';
Generate alignDB from .fas files

    perl init/gen_alignDB.pl -d S288cvsRM11_1a --lt 5000 --parallel 2

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
    [ 'dir_align|da=s', 'dir contains .fas files', { default => $conf->{generate}{dir_align} }, ],
    [ 'length|lt=i',    'threshold of lengths',    { default => $conf->{generate}{length} }, ],
    [ 'parallel=i',     'run in parallel mode',    { default => $conf->{generate}{parallel} }, ],
    [ 'outgroup|o', 'alignments have an outgroup', ],
    { show_defaults => 1, }
    );

$usage->die if $opt->{help};

# record config
$stopwatch->record_conf($opt);

# DBI Data Source Name
my $dsn = sprintf "dbi:mysql:database=%s;host=%s;port=%s", $opt->{db}, $opt->{server}, $opt->{port};

#----------------------------------------------------------#
# Search for all files and push their paths to @files
#----------------------------------------------------------#
$stopwatch->start_message("Generate [$opt->{db}] from directory [$opt->{dir_align}]...");

$opt->{dir_align} = path( $opt->{dir_align} )->stringify;
my @files = sort File::Find::Rule->file->name('*.fas')->in( $opt->{dir_align} );
printf "\n----Total .fas Files: %4s----\n\n", scalar @files;
if ( scalar @files == 0 ) {
    @files = sort File::Find::Rule->file->name('*.fas.gz')->in( $opt->{dir_align} );
    printf "\n----Total .fas.gz Files: %4s----\n\n", scalar @files;
}

#----------------------------------------------------------#
# worker
#----------------------------------------------------------#
my $worker = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;
    my $infile = $chunk_ref->[0];
    my $wid    = MCE->wid;

    my $inner    = AlignDB::Stopwatch->new;
    my $basename = path($infile)->basename;
    $inner->block_message("Process task [$chunk_id] by worker #$wid. [$basename]");

    my $alignDB;
    if ( !$opt->{outgroup} ) {
        $alignDB = AlignDB->new(
            dsn    => $dsn,
            user   => $opt->{username},
            passwd => $opt->{password},
        );
    }
    else {
        $alignDB = AlignDB::Outgroup->new(
            dsn    => $dsn,
            user   => $opt->{username},
            passwd => $opt->{password},
        );
    }

    $alignDB->parse_fas_file( $infile, { threshold => $opt->{length}, } );

    $inner->block_message( "[$basename] has been processed.", "duration" );

    return;
};

#----------------------------------------------------------#
# start insert
#----------------------------------------------------------#
my $mce = MCE->new( max_workers => $opt->{parallel}, chunk_size => 1, );
$mce->foreach( [ sort @files ], $worker );    # foreach also implies chunk_size => 1

$stopwatch->end_message( "All files have been processed.", "duration" );

# store program's meta info to database
AlignDB->new(
    dsn    => $dsn,
    user   => $opt->{username},
    passwd => $opt->{password},
)->add_meta_stopwatch($stopwatch);

exit;

__END__
