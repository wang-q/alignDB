#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Find::Rule;
use File::Spec;

use FindBin;

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
my $working_dir = ".";

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'          => \$help,
    'man'             => \$man,
    'd|working_dir=s' => \$working_dir,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------#
# Run
#----------------------------#
my $str;
if ( $^O eq 'MSWin32' ) {
    $str = "cmd.bat";
}
else {
    $str = "cmd.sh";
}

my @files = sort File::Find::Rule->file->name($str)->in($working_dir);
printf "\n----Total $str Files: %4s----\n\n", scalar @files;

open my $fh, '>>', File::Spec->catfile( $working_dir, "bac_batch_report.txt" );
for my $file (@files) {
    $stopwatch->block_message($file);
    $file = File::Spec->rel2abs($file);
    if ( $^O eq 'MSWin32' ) {
        system "cmd /c $file";
    }
    else {
        system "sh $file";
    }

    print {$fh} $file, "\n";
}
close $fh;

$stopwatch->end_message("All files have been processed.\n");
exit;

__END__

perl bac_batch.pl -d .
