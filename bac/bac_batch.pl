#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use FindBin;
use List::MoreUtils qw(zip);
use Path::Class;
use Template;

use AlignDB::IntSpan;
use AlignDB::Stopwatch;

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

# running options
my $base_dir = '/home/wangq/pork/Bacteria/multi/tree/fas/';
my $run      = "all";

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    'b|base_dir=s' => \$base_dir,
    'r|run=s'      => \$run,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

# prepare to run tasks in @tasks
my @tasks;
if ( $run eq 'all' ) {
    @tasks = ( 1 .. 8 );
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
# init
#----------------------------------------------------------#
my $tt = Template->new;

my $dispatch = {
    1 =>
        qq{ perl $FindBin::Bin/../multi/fasta_malignDB.pl --dir=[% dir %] --db=[% db %] --parallel=4 },
    2 => qq{ perl $FindBin::Bin/../multi/insert_gc_multi.pl --db=[% db %] },
    3 =>
        qq{ perl $FindBin::Bin/../multi/update_multi_gff.pl --db=[% db %] --gff_file=[% gff %] },
    4 =>
        qq{ perl $FindBin::Bin/../stat/mvar_stat_factory.pl --db=[% db %] --output=[% mvar_xls %] },
    5 =>
        qq{ perl $FindBin::Bin/../stat/multi_stat_factory.pl --db=[% db %] --output=[% multi_xls %] },
    6 =>
        qq{ perl $FindBin::Bin/../stat/mgc_stat_factory.pl --db=[% db %] --output=[% mgc_xls %] },
    7 =>
        qq{ perl $FindBin::Bin/../stat/multi_chart_factory.pl -i=[% multi_xls %] },
    8 =>
        qq{ perl $FindBin::Bin/../stat/mgc_chart_factory.pl -i=[% mgc_xls %] },
};

#----------------------------#
# Run
#----------------------------#

my $dir = dir($base_dir);

for my $subdir ( $dir->children ) {
    next unless $subdir->is_dir;
    my ($db_name) = $subdir->dir_list(-1);
    $stopwatch->block_message($db_name);

    my $gff_file;
    while ( my $file = $subdir->next ) {
        next unless -f $file;
        my $filename = $file->absolute->stringify;
        if ( $filename =~ /\.gff/ ) {
            $gff_file = $filename;
            last;
        }
    }

    my $mvar_xls  = "$db_name.mvar.xls";
    my $multi_xls = "$db_name.multi.xls";
    my $mgc_xls   = "$db_name.mgc.xls";

    # use the dispatch template to generate $cmd
    for my $step (@tasks) {
        if ( $step > 6 ) {
            next if $^O ne 'MSWin32';
        }

        my $cmd;
        $tt->process(
            \$dispatch->{$step},
            {   db        => $db_name,
                dir       => $subdir->absolute->stringify,
                gff       => $gff_file,
                mvar_xls  => $mvar_xls,
                multi_xls => $multi_xls,
                mgc_xls   => $mgc_xls,
            },
            \$cmd
        ) or die Template->error;

        $stopwatch->block_message("Processing Step $step");
        $stopwatch->block_message($cmd);
        system $cmd;
        $stopwatch->block_message("Finish Step $step");
    }
}

print "\n";

$stopwatch->end_message;
exit;

__END__

perl bac_batch.pl --base_dir=
