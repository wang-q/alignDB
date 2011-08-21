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

# target, query init values
my $target_taxon_id = $Config->{taxon}{target_taxon_id};
my $target_name     = $Config->{taxon}{target_name};
my $query_taxon_id  = $Config->{taxon}{query_taxon_id};
my $query_name      = $Config->{taxon}{query_name};

my $target = $target_taxon_id . "," . $target_name;    # target sequence
my $query  = $query_taxon_id . "," . $query_name;      # query sequence

# axt
my $axt_dir = $Config->{taxon}{axt_dir};

# running tasks
my $run = "common";

# run in parallel mode
my $parallel = $Config->{generate}{parallel};

# legnth threshold of align
my $axt_threshold = $Config->{generate}{axt_threshold};

# stat parameter
my $sum_threshold = $Config->{stat}{sum_threshold};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'             => \$help,
    'man'                => \$man,
    'server=s'           => \$server,
    'port=i'             => \$port,
    'db=s'               => \$db_name,
    'username=s'         => \$username,
    'password=s'         => \$password,
    'a|axt_dir=s'        => \$axt_dir,
    't|target=s'         => \$target,
    'q|query=s'          => \$query,
    'e|ensembl=s'        => \$ensembl_db,
    'parallel=i'         => \$parallel,
    'at|axt_threshold=i' => \$axt_threshold,
    'st|sum_threshold=i' => \$sum_threshold,
    'run=s'              => \$run,
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
    @tasks = ( 1 .. 3, 10, 21, 30 .. 32, 40, 41 );
}
elsif ( $run eq 'gene' ) {
    @tasks = ( 1 .. 3, 10, 20, 21, 30 .. 33, 40 .. 42 );
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
    1 => "perl $FindBin::Bin/../init/init_alignDB.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name",
    2 => "perl $FindBin::Bin/../init/gen_alignDB.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " -t=\"$target\""
        . " -q=\"$query\""
        . " -a=$axt_dir"
        . " --length=$axt_threshold"
        . " --parallel=$parallel",
    3 => "perl $FindBin::Bin/../init/update_isw_indel_id.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name",
    10 => "perl $FindBin::Bin/../init/insert_gc.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " --parallel=$parallel",
    20 => "perl $FindBin::Bin/../gene/insert_gene.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " -e=$ensembl_db"
        . " --parallel=$parallel",
    21 => "perl $FindBin::Bin/../init/update_sw_cv.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " --parallel=$parallel",
    30 => "perl $FindBin::Bin/../init/update_feature.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " -e=$ensembl_db"
        . " --parallel=$parallel",
    31 => "perl $FindBin::Bin/../init/update_indel_slippage.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name",
    32 => "perl $FindBin::Bin/../init/update_segment.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name",
    33 => "perl $FindBin::Bin/../init/update_snp_dnds.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name",
    40 => "perl $FindBin::Bin/../stat/common_stat_factory.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " -o=$FindBin::Bin/../stat/$db_name.common.xlsx"
        . " -t=$sum_threshold",
    41 => "perl $FindBin::Bin/../stat/gc_stat_factory.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " -o=$FindBin::Bin/../stat/$db_name.gc.xlsx"
        . " -t=$sum_threshold",
    42 => "perl $FindBin::Bin/../stat/gene_stat_factory.pl"
        . " -s=$server"
        . " --port=$port"
        . " -u=$username"
        . " --password=$password"
        . " -d=$db_name"
        . " -o=$FindBin::Bin/../stat/$db_name.gene.xlsx",
    43 => undef,
};

#----------------------------#
# Run
#----------------------------#

# use the dispatch template to generate $cmd
for my $step (@tasks) {
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

perl two_way_batch.pl -d alignDB -e yeast_58 -t "4932,S288C" -q "285006,RM11" -a F:/S288CvsRM11 -at 10000 -st 1000000 --parallel 4 --run common

perl two_way_batch.pl -d HumanvsChimp -e human_58 -t "9606,Human" -q "9598,Chimp" -a /home/wangq/data/UCSC/Human19vsChimp2 -at 10000 -st 10000000 --parallel 6 --run basic
