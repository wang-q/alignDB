#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use Win32::OLE qw(in);
use Win32::OLE::Const;
use Win32::OLE::Variant;
use Win32::OLE::NLS qw(:LOCALE :DATE);
use Win32::OLE::Const 'Microsoft Excel';

use File::Find::Rule;
use File::Spec;

$Win32::OLE::Warn = 3;    # die on errors...

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
# running options
my $dir = '.';

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'    => \$help,
    'man'       => \$man,
    'd|dir=s' => \$dir,
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
$excel->{DisplayAlerts} = 0;

my @xlsx_files = File::Find::Rule->file->name('*.xlsx')->in($dir);

for my $xlsx (@xlsx_files) {
    $xlsx = File::Spec->rel2abs($xlsx);
    print $xlsx, "\n";
    my $newbook = $excel->Workbooks->Open($xlsx)
        or die "Cannot open xlsx file\n";
    my $xls = $xlsx;
    $xls =~ s/xlsx$/xls/i;
    
    # xlExcel8    56 Excel8 
    # xlExcel9795 43 Excel9795 # this doesn't work
    $newbook->SaveAs($xls, 56);
}

$excel->Quit;
__END__
