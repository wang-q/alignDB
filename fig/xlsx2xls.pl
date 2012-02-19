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

my $suffix = "*.xlsx";

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'    => \$help,
    'man'       => \$man,
    'd|dir=s' => \$dir,
    's|suffix=s' => \$suffix,
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

my @xlsx_files = File::Find::Rule->file->name($suffix)->in($dir);

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

=head1 NAME

    xlsx2xls.pl - convert xlsx to xls

=head1 SYNOPSIS
    perl xlsx2xls.pl -d stat -s *.chart.xlsx

    xlsx2xls.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        -d, --dir           dir containing xlsx files
        -s, --suffix        suffix, default is *.xlsx

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do someting
useful with the contents thereof.

=cut

