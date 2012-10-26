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

use FindBin;
use lib "$FindBin::Bin/../lib";
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

# run in parallel mode
my $parallel = 1;

my $axt_threshold = $Config->{generate}{axt_threshold};
my $sum_threshold = $Config->{stat}{sum_threshold};

# running options
my $bz = "$FindBin::Bin/../../blastz/bz.pl";
my $pair_file;

my $dir_as_taxon;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'             => \$help,
    'man'                => \$man,
    'b|bz=s'             => \$bz,
    'p|parallel=i'       => \$parallel,
    'at|axt_threshold=i' => \$axt_threshold,
    'st|sum_threshold=i' => \$sum_threshold,
    'f|pair_file=s'      => \$pair_file,
    'd|dir_as_taxon=s'   => \$dir_as_taxon,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
my ( $volume, $directories, undef ) = File::Spec->splitpath($pair_file);
my $working_dir = File::Spec->catpath( $volume, $directories );
my $report = File::Spec->catfile( $working_dir, "pair_db.txt" );
open my $fh, '>', $report;

my $tt = Template->new;

my $dispatch = {
    0 => "perl $bz"
        . " -dt [% tfile %] -dq [% qfile %] -dl [% ldir %]"
        . " -s set01 --parallel [% parallel %] -pb lastz --lastz",
    1 => "perl $FindBin::Bin/../init/init_alignDB.pl" . " --db [% db %] ",
    2 => "perl $FindBin::Bin/../init/gen_alignDB.pl"
        . " --db [% db %]"
        . " -a [% ldir %] [% tq %]"
        . " --length [% at %]"
        . " --parallel [% parallel %]",
    21 => "perl $FindBin::Bin/../init/update_sw_cv.pl"
        . " --db [% db %]"
        . " --parallel [% parallel %]",
    40 => "perl $FindBin::Bin/../stat/common_stat_factory.pl"
        . " --db [% db %]"
        . " -o [% common_file %]"
        . " -threshold [% st %]",
};

#----------------------------#
# Run
#----------------------------#

my @lines = read_file($pair_file);
for (@lines) {
    chomp;
    my ( $tfile, $qfile ) = split /,/, $_;
    $tfile or next;

    my $t_base = basename($tfile);
    $t_base =~ s/\..+?$//;
    my $q_base = basename($qfile);
    $q_base =~ s/\..+?$//;
    my $db = "${t_base}vs${q_base}";
    my $ldir = File::Spec->catdir( $working_dir, $db );
    print {$fh} "$db\n";

    my $tq;
    if ($dir_as_taxon) {
        $tq = " -t=\"$t_base,$t_base\"" . " -q=\"$q_base,$q_base\"";
    }

    my $common_file = File::Spec->catfile( $working_dir, "$db.common.xlsx" );

    # use the dispatch template to generate $cmd
    for my $step ( 0 .. 2, 21, 40 ) {

        my $cmd;
        $tt->process(
            \$dispatch->{$step},
            {   tfile       => $tfile,
                qfile       => $qfile,
                ldir        => $ldir,
                db          => $db,
                tq          => $tq,
                parallel    => $parallel,
                common_file => $common_file,
                at          => $axt_threshold,
                st          => $sum_threshold,
            },
            \$cmd
        ) or die Template->error;

        $stopwatch->block_message("Processing Step $step");
        $stopwatch->block_message($cmd);
        system $cmd;
        $stopwatch->block_message("Finish Step $step");
    }
}
close $fh;
print "\n";

$stopwatch->end_message;
exit;

__END__

perl seq_pair_batch.pl -d 1 -p 2 -f "e:\wq\Scripts\alignDB\bac\seq_pair.csv" 
