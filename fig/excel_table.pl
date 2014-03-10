#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use Path::Class;
use Set::Scalar;

use Win32::OLE qw(in);
use Win32::OLE::Const;
use Win32::OLE::Variant;
use Win32::OLE::NLS qw(:LOCALE :DATE);
use Win32::OLE::Const 'Microsoft Excel';

$Win32::OLE::Warn = 2;    # die on errors...

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
# running options
my $file_yaml = 'Fig.S1.yaml';

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'    => \$help,
    'man'       => \$man,
    'i|input=s' => \$file_yaml,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# init
#----------------------------------------------------------#
# The structure of yaml file:
# filename:sheetname:chart_serial:[x, y] coordinates
my $dispatch = LoadFile($file_yaml);
my $ranges   = $dispatch->{ranges};
my $texts    = $dispatch->{texts};
my $borders  = $dispatch->{borders};

# Excel files should be located in the same dir as the yaml file.
# new excel table file will be named after yaml file.
$file_yaml = file($file_yaml)->absolute;
my $base_dir = $file_yaml->dir->stringify;
my $newfile  = $file_yaml->stringify;
$newfile =~ s/\.ya?ml$/.xlsx/;
unlink $newfile if -e $newfile;

my $excel;    # excel object
unless ( $excel = Win32::OLE->new("Excel.Application") ) {
    die "Cannot init Excel.Application\n";
}

my $newbook = $excel->Workbooks->Add
    or die "Cannot create xls file\n";
$newbook->SaveAs($newfile);

for ( in $newbook->Worksheets ) {
    $_->Delete;
    $newbook->Worksheets->Count < 2 and last;
}

my $sheetname = $file_yaml->basename;
$sheetname =~ s/[^\w]/_/g;
if ( length $sheetname > 30 ) {
    $sheetname = substr $sheetname, 0, 30;
}
my $newsheet = $newbook->Worksheets(1);
$newsheet->{Name} = $sheetname;

#----------------------------------------------------------#
# Paste every ranges
#----------------------------------------------------------#
for my $filename ( sort keys %{$ranges} ) {
    printf "[file: %s]\n", $filename;

    # open xls file
    my $efile = file( $base_dir, $filename )->stringify;
    if ( !-e $efile ) {
        warn "File not exists: $efile\n";
        next;
    }
    my $workbook;
    unless ( $workbook = $excel->Workbooks->Open($efile) ) {
        die "Cannot open xls file\n";
    }
    for my $sheetname ( sort keys %{ $ranges->{$filename} } ) {
        printf "[sheet: %s]\n", $sheetname;

        # get sheet
        my @sheets;
        for my $sheet ( in $workbook->Worksheets ) {
            push @sheets, $sheet->{Name};
        }
        my $name_set = Set::Scalar->new(@sheets);

        if ( !$name_set->has($sheetname) ) {
            print " " x 4, "sheet not exists!\n";
            next;
        }
        my $sheet = $workbook->Worksheets($sheetname);
        for my $range ( @{ $ranges->{$filename}{$sheetname} } ) {
            printf "[range]\n";

            # copy
            my $copy_range = $sheet->Range( $range->{copy} );
            $copy_range->Copy;

            # paste
            my $paste_range = $newsheet->Range( $range->{paste} );
            $paste_range->Cells(1)->PasteSpecial;
            $paste_range->Merge if $range->{merge};

            # clear clipboard
            $excel->{CutCopyMode} = 0;
        }
    }
    $workbook->Close;
}

#----------------------------------------------------------#
# Write every texts
#----------------------------------------------------------#
for my $i ( @{$texts} ) {
    my $text = $i->{text};
    my $pos  = $i->{pos};

    my $range = $newsheet->Range($pos);
    $range->Cells(1)->{Value} = $text;

    $range->Cells(1)->{Font}->{Name} = $i->{font} || "Arial";
    $range->Cells(1)->{Font}->{Size} = $i->{font} || 10;
    $range->Cells(1)->{Font}->{Bold}   = $i->{bold}   ? 1 : 0;
    $range->Cells(1)->{Font}->{Italic} = $i->{italic} ? 1 : 0;

    $range->Merge if $i->{merge};

    printf "[text: %s]\n", $text;
}

#----------------------------------------------------------#
# Borders
#----------------------------------------------------------#
for my $i ( @{$borders} ) {
    my $range = $i->{range};

    $range = $newsheet->Range($range);
    if ( $i->{top} ) {
        $range->Borders(xlEdgeTop)->{LineStyle} = xlContinuous;
    }
    if ( $i->{bottom} ) {
        $range->Borders(xlEdgeBottom)->{LineStyle} = xlContinuous;
    }
    if ( $i->{left} ) {
        $range->Borders(xlEdgeLeft)->{LineStyle} = xlContinuous;
    }
    if ( $i->{right} ) {
        $range->Borders(xlEdgeRight)->{LineStyle} = xlContinuous;
    }

    printf "[border]\n";
}

#----------------------------#
# Style
#----------------------------#
if ( $dispatch->{autofit} ) {
    $newsheet->{UsedRange}->{Columns}->AutoFit;
}

#----------------------------------------------------------#
# Cleanup and Save
#----------------------------------------------------------#
$newbook->Save;
$excel->Quit;

__END__

=head1 SYNOPSIS

perl excel_table.pl -i Fig.S1.yaml
