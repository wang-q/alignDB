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

# alignments have an outgroup
my $outgroup;

my @gff_files;
my @rm_gff_files;

# alignment
my $dir_align        = $Config->{taxon}{dir_align};
my $length_threshold = $Config->{generate}{length_threshold};

my $block;         # input is galaxy style blocked fasta
my $file_id_of;    # taxon_id-name mapping file

# running tasks
my $run = "common";

# run in parallel mode
my $parallel = $Config->{generate}{parallel};

# number of alignments process in one child process
my $batch_number = $Config->{generate}{batch};

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'                 => \$help,
    'man'                    => \$man,
    's|server=s'             => \$server,
    'P|port=i'               => \$port,
    'u|username=s'           => \$username,
    'p|password=s'           => \$password,
    'd|db=s'                 => \$db_name,
    'o|outgroup'             => \$outgroup,
    'da|dir_align=s'         => \$dir_align,
    'e|ensembl=s'            => \$ensembl_db,
    'gff_files=s'            => \@gff_files,
    'rm_gff_files=s'         => \@rm_gff_files,
    'id|id_of=s'             => \$file_id_of,
    'block'                  => \$block,
    'parallel=i'             => \$parallel,
    'batch=i'                => \$batch_number,
    'lt|length_threshold=i'  => \$length_threshold,
    'r|run=s'                => \$run,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

# prepare to run tasks in @tasks
my @tasks;
if ( $run eq 'all' ) {
    @tasks = ( 1 .. 45 );
}
elsif ( $run eq 'skeleton' ) {
    @tasks = ( 1 .. 3 );
}
elsif ( $run eq 'basic' ) {
    @tasks = ( 1 .. 5, 40 );
}
elsif ( $run eq 'common' ) {
    @tasks = ( 1 .. 5, 21, 30, 31, 40, 41, 44 );
}
elsif ( $run eq 'gc' ) {
    @tasks = ( 1 .. 5, 10, 21, 30 .. 32, 40 .. 42 );
}
elsif ( $run eq 'gene' ) {
    @tasks = ( 1 .. 5, 10, 20, 21, 30 .. 33, 40 .. 44 );
}
elsif ( $run eq 'stat' ) {
    @tasks = ( 40 .. 44 );
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

my $gff_file = join ",", @gff_files;
my $rm_gff_file = join ",", @rm_gff_files;

#----------------------------------------------------------#
# dispatch table
#----------------------------------------------------------#
my $dispatch = {
    1 => "perl $FindBin::Bin/../init/init_alignDB.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name",
    2 => "perl $FindBin::Bin/../init/gen_alignDB_fas.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " --db $db_name"
        . " --da $dir_align"
        . " -lt $length_threshold"
        . " --parallel $parallel"
        . " --batch $batch_number"
        . " --id $file_id_of"
        . ( $outgroup ? " --outgroup" : "" )
        . ( $block    ? " --block"    : "" ),
    3 => undef,
    5 => "perl $FindBin::Bin/../init/insert_isw.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name"
        . " --parallel $parallel"
        . " --batch $batch_number"
        . ( $outgroup ? " --outgroup" : "" ),
    10 => "perl $FindBin::Bin/../init/insert_gc.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name"
        . " --parallel $parallel"
        . " --batch $batch_number",
    20 => "perl $FindBin::Bin/../gene/insert_gene.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name"
        . " -e $ensembl_db"
        . " --parallel $parallel",
    21 => "perl $FindBin::Bin/../init/update_sw_cv.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name"
        . " --parallel $parallel"
        . " --batch $batch_number",
    30 => "perl $FindBin::Bin/../init/update_feature.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name"
        . " -e $ensembl_db"
        . " --parallel $parallel"
        . " --batch $batch_number",
    '30gff' => "perl $FindBin::Bin/../init/update_feature_gff.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name"
        . " --parallel $parallel"
        . " --batch $batch_number"
        . " --gff_files $gff_file"
        . " --rm_gff_files $rm_gff_file",
    31 => "perl $FindBin::Bin/../init/update_indel_slippage.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name",
    32 => "perl $FindBin::Bin/../init/update_segment.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name",
    33 => "perl $FindBin::Bin/../gene/update_snp_dnds.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name",
    40 => "perl $FindBin::Bin/../stat/common_stat_factory.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name"
        . " -o $db_name.common.xlsx",
    41 => "perl $FindBin::Bin/../stat/multi_stat_factory.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name"
        . " -o $db_name.multi.xlsx",
    42 => "perl $FindBin::Bin/../stat/gc_stat_factory.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name"
        . " -o $db_name.gc.xlsx",
    43 => "perl $FindBin::Bin/../stat/gene_stat_factory.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name"
        . " -o $db_name.gene.xlsx",
    44 => "perl $FindBin::Bin/../stat/mvar_stat_factory.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name"
        . " -o $db_name.mvar.xlsx",
};

#----------------------------#
# Run
#----------------------------#

# use the dispatch template to generate $cmd
for my $step (@tasks) {
    if ( @gff_files and $step == 30 ) {
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

=head1 SYNOPSIS

perl multi_way_batch.pl -d S288CvsThree_10k -e yeast_58 -da F:/S288CvsThree_10k -lt 5000 -st 0 --parallel=6 --run all

perl multi_way_batch.pl -d S288CvsTen_10k -e yeast_58 -da F:/S288CvsTen_10k -lt 5000 -st 0 --parallel=6 --run all

perl multi_way_batch.pl -d S288CvsSix_10k -r stat

perl multi_way_batch.pl -d AthvsFive -da d:\data\alignment\arabidopsis\AthvsFive\ -lt 5000 -st 0 --parallel 4 --run 1-3,21,40

perl ~/Scripts/alignDB/extra/multi_way_batch.pl -d AthvsV_mafft --block --id 3702 -da ~/data/alignment/arabidopsis19/AthvsV_mafft -lt 5000 -st 0 --parallel 4 --run 1-3,21,40

