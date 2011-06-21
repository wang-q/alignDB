#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

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
my $server     = $Config->{database}{server};
my $port       = $Config->{database}{port};
my $username   = $Config->{database}{username};
my $password   = $Config->{database}{password};
my $db_name    = $Config->{database}{db};
my $ensembl_db = $Config->{database}{ensembl};

my $gff_file = '';

# axt
my $fasta_dir        = '/home/wangq/Date/Alignment/yeast6/';
my $length_thredhold = $Config->{ref}{length_threshold};

# running tasks
my $run = "all";

# run in parallel mode
my $parallel = $Config->{generate}{parallel};

my $all_freq = 3;

# stat parameter
my $sum_threshold = $Config->{stat}{sum_threshold};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'                => \$help,
    'man'                   => \$man,
    'server=s'              => \$server,
    'port=i'                => \$port,
    'db=s'                  => \$db_name,
    'username=s'            => \$username,
    'password=s'            => \$password,
    'f|fasta_dir=s'         => \$fasta_dir,
    'e|ensembl=s'           => \$ensembl_db,
    'gff_file=s'            => \$gff_file,
    'parallel=i'            => \$parallel,
    'lt|length_thredhold=i' => \$length_thredhold,
    'st|sum_threshold=i'    => \$sum_threshold,
    'all_freq=s'            => \$all_freq,
    'run=s'                 => \$run,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

# prepare to run tasks in @tasks
my @tasks;
if ( $run eq 'all' ) {
    @tasks = ( 1 .. 43 );
}
elsif ( $run eq 'basic' ) {
    @tasks = ( 1 .. 3 );
}
elsif ( $run eq 'common' ) {
    @tasks = ( 1 .. 3, 30, 31, 40 );
}
elsif ( $run eq 'gc' ) {
    @tasks = ( 1 .. 3, 10, 11, 30 .. 32, 40, 41 );
}
elsif ( $run eq 'gene' ) {
    @tasks = ( 1 .. 3, 10, 11, 20, 21, 30 .. 33, 40 .. 42 );
}
elsif ( $run eq 'stat' ) {
    @tasks = ( 40 .. 42 );
}
else {
    $run =~ s/\"\'//s;
    if ( AlignDB::IntSpan->valid($run) ) {
        my $set = AlignDB::IntSpan->new($run);
        @tasks = $set->elements;
    }
    else {
        @tasks = grep {/\d/} split /\s/, $run;
    }
}

#----------------------------------------------------------#
# dispatch table
#----------------------------------------------------------#
my $dispatch = {
    1 => "perl $FindBin::Bin/../multi/fasta_malignDB.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " --fasta_dir=$fasta_dir"
        . " --length=$length_thredhold"
        . " --parallel=$parallel",
    2  => undef,
    3  => undef,
    10 => "perl $FindBin::Bin/../multi/insert_gc_multi.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " --parallel=$parallel",
    11 => "perl $FindBin::Bin/../init/update_isw_cv.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name",
    20 => "perl $FindBin::Bin/../gene/insert_gene.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " -e=$ensembl_db"
        . " --parallel=$parallel"
        . " --multi",
    21 => "perl $FindBin::Bin/../gene/update_sw_cv.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " --parallel=$parallel"
        . " --multi",
    30 => "perl $FindBin::Bin/../multi/update_multi.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " -e=$ensembl_db",
    '30gff' => "perl $FindBin::Bin/../multi/update_multi_gff.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " --gff_file=$gff_file",
    31 => undef,
    32 => undef,
    33 => "perl $FindBin::Bin/../init/update_snp_dnds.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " --multi",
    40 => "perl $FindBin::Bin/../stat/multi_stat_factory.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " -o=$FindBin::Bin/../stat/$db_name.multi.xlsx"
        . " --freq=$all_freq",
    41 => "perl $FindBin::Bin/../stat/mgc_stat_factory.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " -o=$FindBin::Bin/../stat/$db_name.mgc.xlsx"
        . " -t=$sum_threshold",
    42 => "perl $FindBin::Bin/../stat/gene_stat_factory.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " -o=$FindBin::Bin/../stat/$db_name.gene.xlsx",
    43 => "perl $FindBin::Bin/../stat/mvar_stat_factory.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " -o=$FindBin::Bin/../stat/$db_name.mvar.xlsx",
};

#----------------------------#
# Run
#----------------------------#

# use the dispatch template to generate $cmd
for my $step (@tasks) {
    if ( $gff_file and $step == 30 ) {
        $step = '30gff';
    }

    my $cmd = $dispatch->{$step};

    next unless $cmd;

    $stopwatch->block_message("Processing Step $step");
    $stopwatch->block_message($cmd);
    system $cmd;
    $stopwatch->block_message("Finish Step $step");
}

print "\n";

$stopwatch->end_message;
exit;

__END__

perl multi_way_batch.pl -d S288CvsThree_10k -e yeast_58 -f F:/S288CvsThree_10k --all_freq 3 -lt 10000 -st 100000 --parallel=6 --run all

perl multi_way_batch.pl -d S288CvsTen_10k -e yeast_58 -f F:/S288CvsTen_10k --all_freq 10 -lt 10000 -st 100000 --parallel=6 --run all
