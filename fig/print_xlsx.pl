#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use Win32::OLE qw(in with);
use Win32::OLE::Const;
use Win32::OLE::Variant;
use Win32::OLE::NLS qw(:LOCALE :DATE);
use Win32::OLE::Const 'Microsoft Excel';

use Path::Class;

$Win32::OLE::Warn = 2;    # die on errors...

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
# running options
my $file_xlsx = 'Fig.S1.yaml';

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'    => \$help,
    'man'       => \$man,
    'i|input=s' => \$file_xlsx,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
my $excel;    # excel object
unless ( $excel = Win32::OLE->new("Excel.Application") ) {
    die "Cannot init Excel.Application\n";
}

my $workbook  = $excel->Workbooks->Open($file_xlsx);
my $sheet     = $workbook->Worksheets(1);
with(
    $sheet->PageSetup,
    Orientation    => xlLandscape,
    FitToPagesWide => 1,
    FitToPagesTall => 1,
    PaperSize      => xlPaperA4,
    Zoom           => Variant(VT_BOOL, 0),
);
my $pobj = $excel->ActivePrinter;
if ($pobj) {
    print "start - Name of default printer => $pobj \n";
}
$pobj = $excel->ActivePrinter;
if ($pobj) {
    print "after setting - Name of default printer => $pobj \n";
}
$workbook->PrintOut();
$workbook->Close(0);

$excel->Quit;
