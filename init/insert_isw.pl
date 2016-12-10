#!/usr/bin/env perl
use strict;
use warnings;
use autodie;

use Getopt::Long::Descriptive;
use Config::Tiny;
use FindBin;
use YAML::Syck;

use MCE;

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
Update indel-sliding windows

    perl init/insert_isw.pl -d S288cvsRM11_1a --parallel 2

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
    [ 'parallel=i', 'run in parallel mode',       { default => $conf->{generate}{parallel} }, ],
    [ 'batch=i',    '#alignments in one process', { default => $conf->{generate}{batch} }, ],
    [ 'outgroup|o', 'alignments have an outgroup', ],
    { show_defaults => 1, }
    );

$usage->die if $opt->{help};

# record config
$stopwatch->record_conf($opt);

# DBI Data Source Name
my $dsn = sprintf "dbi:mysql:database=%s;host=%s;port=%s", $opt->{db}, $opt->{server}, $opt->{port};

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
$stopwatch->start_message("Update isw-indel relationship of [$opt->{db}]...");

my @jobs;
{
    my $alignDB = AlignDB->new(
        dsn    => $dsn,
        user   => $opt->{username},
        passwd => $opt->{password},
    );

    # empty tables
    print "Emptying tables...\n";
    $alignDB->empty_table( 'isw', );

    @jobs = @{ $alignDB->get_align_ids };
}

#----------------------------------------------------------#
# worker
#----------------------------------------------------------#
my $worker = sub {
    my ( $self, $chunk_ref, $chunk_id ) = @_;
    my @align_ids = @{$chunk_ref};
    my $wid       = MCE->wid;

    $stopwatch->block_message("Process task [$chunk_id] by worker #$wid");

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

    # for each alignment
    for my $align_id (@align_ids) {
        $alignDB->process_message($align_id);
        $alignDB->insert_isw($align_id);
        $alignDB->isw_snp_fk($align_id);
        $alignDB->update_D_values($align_id) if $opt->{outgroup};
    }

    return;
};

#----------------------------------------------------------#
# start update
#----------------------------------------------------------#
my $mce = MCE->new( max_workers => $opt->{parallel}, chunk_size => $opt->{batch}, );
$mce->forchunk( \@jobs, $worker, );

$stopwatch->end_message;

# store program's meta info to database
AlignDB->new(
    dsn    => $dsn,
    user   => $opt->{username},
    passwd => $opt->{password},
)->add_meta_stopwatch($stopwatch);

exit;

__END__
