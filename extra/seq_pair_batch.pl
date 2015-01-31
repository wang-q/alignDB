#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Basename;
use File::Slurp;
use File::Spec;
use Template;

use AlignDB::Run;
use AlignDB::IntSpan;
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

my $length_threshold = $Config->{generate}{length_threshold};

# executable file location
my $kent_bin = "~/bin/x86_64";

# running options
my $bz_path = "$FindBin::Bin/../../egaz";
my $pair_file;

my $dir_as_taxon;    # use dir name as taxon_id
my $file_id_of;      # taxon_id-name mapping file

# running tasks
my $task = "0-2,21,40";

# run in parallel mode
my $parallel = $Config->{generate}{parallel};

# number of alignments process in one child process
my $batch_number = $Config->{generate}{batch};

my $init_taxon = "$FindBin::Bin/../data/taxon.csv";
my $init_chr   = "$FindBin::Bin/../data/chr_length.csv";

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'                => \$help,
    'man'                   => \$man,
    's|server=s'            => \$server,
    'P|port=i'              => \$port,
    'u|username=s'          => \$username,
    'p|password=s'          => \$password,
    'bz=s'                  => \$bz_path,
    'bin|kent_bin=s'        => \$kent_bin,
    'parallel=i'            => \$parallel,
    'batch=i'               => \$batch_number,
    'lt|length_threshold=i' => \$length_threshold,
    'f|pair_file=s'         => \$pair_file,
    'dir_as_taxon'          => \$dir_as_taxon,
    'id|id_of=s'            => \$file_id_of,
    'r|run=s'               => \$task,
    'taxon|init_taxon=s'    => \$init_taxon,
    'chr|init_chr=s'        => \$init_chr,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

my @tasks;
{
    $task =~ s/\"\'//s;
    if ( AlignDB::IntSpan->valid($task) ) {
        my $set = AlignDB::IntSpan->new($task);
        @tasks = $set->elements;
    }
    else {
        @tasks = grep {/\d/} split /\s/, $task;
    }
}

# The server has 16 cpu cores
if ($parallel > 4) {
    $parallel = 4;
}

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
my ( $volume, $directories, undef ) = File::Spec->splitpath($pair_file);
my $working_dir = File::Spec->catpath( $volume, $directories );

my @lines = read_file($pair_file);
my $id_of = {};
if ($file_id_of) {
    open my $fh, '<', $file_id_of;
    while (<$fh>) {
        chomp;
        my ( $id, $name ) = split /,/;
        $id_of->{$name} = $id;
    }
    close $fh;
}

#----------------------------#
# worker
#----------------------------#
my $worker = sub {
    my $line = shift;
    chomp $line;
    my ( $tfile, $qfile ) = split /,/, $line;
    $tfile or return;

    my $t_base = basename($tfile);
    $t_base =~ s/\..+?$//;
    my $q_base = basename($qfile);
    $q_base =~ s/\..+?$//;
    my $db = "${t_base}vs${q_base}";
    my $ldir = File::Spec->catdir( $working_dir, $db );

    my $tq;
    if ($dir_as_taxon) {
        $tq = sprintf " -t=\"%s,%s\" -q=\"%s,%s\"", $t_base, $t_base, $q_base,
            $q_base;
    }
    elsif ($file_id_of) {
        $tq = sprintf " -t=\"%s,%s\" -q=\"%s,%s\"", $id_of->{$t_base}, $t_base,
            $id_of->{$q_base}, $q_base;
    }
    else {
        warn "Taxon info not specified.\n";
    }

    my $common_file = File::Spec->catfile( $working_dir, "$db.common.xlsx" );

    my $tt = Template->new;

    my $dispatch = {
        0 => "perl $bz_path/bz.pl"
            . " -dt [% tfile %] -dq [% qfile %] -dl [% ldir %]"
            . " -s set01 --parallel [% parallel %] -pb lastz --lastz",
        1 => "perl $FindBin::Bin/../init/init_alignDB.pl"
            . " -s $server"
            . " --port $port"
            . " -u $username"
            . " --password $password"
            . " -d $db"
            . ( $init_taxon ? " -taxon $init_taxon" : "" )
            . ( $init_chr   ? " -chr $init_chr"     : "" ),
        2 => "perl $FindBin::Bin/../init/gen_alignDB.pl"
            . " -s $server"
            . " --port $port"
            . " -u $username"
            . " --password $password"
            . " -d $db"
            . " -da [% ldir %]"
            . " [% tq %]"
            . " --length   [% lt %]"
            . " --parallel $parallel",
        5 => "perl $FindBin::Bin/../init/insert_isw.pl"
            . " -s $server"
            . " --port $port"
            . " -u $username"
            . " --password $password"
            . " -d $db"
            . " --parallel $parallel"
            . " --batch $batch_number",
        21 => "perl $FindBin::Bin/../init/update_sw_cv.pl"
            . " -s $server"
            . " --port $port"
            . " -u $username"
            . " --password $password"
            . " -d $db"
            . " --parallel $parallel"
            . " --batch $batch_number",
        40 => "perl $FindBin::Bin/../stat/common_stat_factory.pl"
            . " -s $server"
            . " --port $port"
            . " -u $username"
            . " --password $password"
            . " -d $db"
            . " -o [% common_file %]",
        100 => "perl $bz_path/bz.pl"
            . " -dt [% tfile %] -dq [% qfile %] -dl [% ldir %]"
            . " -s set01 --parallel [% parallel %] --noaxt -pb lastz --lastz",
        101 => "perl $bz_path/lpcna.pl"
            . " -bin $kent_bin"
            . " -dt [% tfile %] -dq [% qfile %] -dl [% ldir %]"
            . " --parallel [% parallel %]",
        102 => "perl $bz_path/amp.pl"
            . " -bin $kent_bin"
            . " -dt [% tfile %] -dq [% qfile %] -dl [% ldir %]"
            . " -syn --parallel [% parallel %]",
    };

    # use the dispatch template to generate $cmd
    for my $step (@tasks) {

        my $cmd;
        $tt->process(
            \$dispatch->{$step},
            {   tfile       => $tfile,
                qfile       => $qfile,
                ldir        => $ldir,
                db          => $db,
                tq          => $tq,
                parallel    => $parallel,
                batch       => $batch_number,
                common_file => $common_file,
                lt          => $length_threshold,
            },
            \$cmd
        ) or die Template->error;

        $stopwatch->block_message("Processing Step $step");
        $stopwatch->block_message($cmd);
        system $cmd;
        $stopwatch->block_message("Finish Step $step");
    }
    print "\n";
};

#----------------------------------------------------------#
# start
#----------------------------------------------------------#
my $run = AlignDB::Run->new(
    parallel => $parallel,
    jobs     => \@lines,
    code     => $worker,
);
$run->run;

$stopwatch->end_message("All files have been processed.");
exit;

__END__

perl seq_pair_batch.pl -d 1 -p 2 -f "e:\wq\Scripts\alignDB\bac\seq_pair.csv" 
