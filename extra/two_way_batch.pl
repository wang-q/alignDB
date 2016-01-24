#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long qw(HelpMessage);
use Config::Tiny;
use FindBin;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

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

two_way_batch.pl - Batch process two-way alignDB

=cut

GetOptions(
    'help|?' => sub { HelpMessage(0) },
    'server|s=s'            => \( my $server           = $Config->{database}{server} ),
    'port|P=i'              => \( my $port             = $Config->{database}{port} ),
    'db|d=s'                => \( my $db_name          = $Config->{database}{db} ),
    'username|u=s'          => \( my $username         = $Config->{database}{username} ),
    'password|p=s'          => \( my $password         = $Config->{database}{password} ),
    'ensembl|e=s'           => \( my $ensembl_db       = $Config->{database}{ensembl} ),
    'dir_align|da=s'        => \( my $dir_align        = $Config->{taxon}{dir_align} ),
    'target|t=s'            => \( my $target_name      = $Config->{taxon}{target_name} ),
    'query|q=s'             => \( my $query_name       = $Config->{taxon}{query_name} ),
    'gff_files=s'           => \my @gff_files,
    'rm_gff_files=s'        => \my @rm_gff_files,
    'parallel=i'            => \( my $parallel         = $Config->{generate}{parallel} ),
    'batch=i'               => \( my $batch_number     = $Config->{generate}{batch} ),
    'length_threshold|lt=i' => \( my $length_threshold = $Config->{generate}{length_threshold} ),
    'run|r=s' => \( my $run        = "common" ),                                     # running tasks
    'chr=s'   => \( my $init_chr   = "$FindBin::RealBin/../data/chr_length.csv" ),
) or HelpMessage(1);

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

my $gff_file    = join ",", @gff_files;
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
        . " -d $db_name"
        . ( $init_chr   ? " --chr $init_chr"     : "" ),
    2 => "perl $FindBin::Bin/../init/gen_alignDB.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " --db $db_name"
        . " -t $target_name"
        . " -q $query_name"
        . " --da $dir_align"
        . " -lt $length_threshold"
        . " --parallel $parallel",
    3 => undef,
    5 => "perl $FindBin::Bin/../init/insert_isw.pl"
        . " -s $server"
        . " --port $port"
        . " -u $username"
        . " --password $password"
        . " -d $db_name"
        . " --parallel $parallel"
        . " --batch $batch_number",
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
    41 => undef,
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
