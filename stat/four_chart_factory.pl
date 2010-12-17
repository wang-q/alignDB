#!/usr/bin/perl
use warnings;
use strict;

use Getopt::Long;
use Pod::Usage;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB::Excel;
use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $infile        = '';
my $outfile       = '';
my $jc_correction = 0;
my $time_stamp    = 1;
my $add_index_sheet    = 1;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    'infile=s'     => \$infile,
    'outfile=s'    => \$outfile,
    'jc=s'         => \$jc_correction,
    'time_stamp=s' => \$time_stamp,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init section
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new();
$stopwatch->start_message("Processing $infile...");

my $excel_obj;
if ($outfile) {
    $excel_obj = AlignDB::Excel->new(
        infile  => $infile,
        outfile => $outfile,
    );
}
else {
    $excel_obj = AlignDB::Excel->new( infile => $infile, );
    $outfile = $excel_obj->outfile;
}

#----------------------------------------------------------#
# START
#----------------------------------------------------------#
# jc
$excel_obj->jc_correction if $jc_correction;

#----------------------------------------------------------#
# draw charts section
#----------------------------------------------------------#
{

    #----------------------------#
    # worksheet -- Dxr
    #----------------------------#
    my @sheets = qw{
        indel_onc_nonslip_dxr
        indel_oc_nonslip_dxr
        indel_onc_nonslip_dxr_1-50
        indel_oc_nonslip_dxr_1-50
    };

    foreach (@sheets) {
        my $sheet_name = $_;
        my %option     = (
            chart_serial  => 1,
            x_column      => 1,
            y_column      => 2,
            y_last_column => 2,
            first_row     => 2,
            last_row      => 32,
            x_max_scale   => 30,
            y_scale_unit  => 0.001,
            x_title       => "Distance to indel (d1)",
            y_title       => "Nucleotide diversity",
            Height        => 283.7,
            Width         => 453.9,
            Top           => 12.75,
            Left          => 360,
        );

        $excel_obj->draw_y( $sheet_name, \%option );

        # chart 2
        $option{chart_serial}++;
        $option{y_column}      = 3;
        $option{y_last_column} = 3;
        $option{Top} += $option{Height} + 12.75;

        $excel_obj->draw_y( $sheet_name, \%option );

        # chart 3
        $option{chart_serial}++;
        $option{y_column}      = 4;
        $option{y_last_column} = 4;
        $option{Top} += $option{Height} + 12.75;

        $excel_obj->draw_y( $sheet_name, \%option );

        # chart 4
        $option{chart_serial}++;
        $option{y_column}      = 5;
        $option{y_last_column} = 5;
        $option{Top} += $option{Height} + 12.75;

        $excel_obj->draw_y( $sheet_name, \%option );
    }

}

{

    #----------------------------#
    # worksheet -- indel_occured, indel_size_didn
    #----------------------------#
    my @sheets = qw{
        indel_onc_distance indel_oc_distance
        indel_onc_nonslip indel_onc_slip
        indel_oc_nonslip indel_oc_slip
        indel_onc_nonslip_1-50
        indel_oc_nonslip_1-50
        indel_size_didn_1-50
    };

    foreach (@sheets) {
        my $sheet_name = $_;
        my %option     = (
            chart_serial  => 1,
            x_column      => 1,
            y_column      => 3,
            y_last_column => 4,
            first_row     => 2,
            last_row      => 7,
            x_max_scale   => 5,
            y_scale_unit  => 0.001,
            x_title       => "Distance to indel (d1)",
            y_title       => "Nucleotide diversity",
            Height        => 283.7,
            Width         => 453.9,
            Top           => 12.75,
            Left          => 440,
        );

        $excel_obj->draw_y( $sheet_name, \%option );

        # chart 2
        $option{chart_serial}++;
        $option{y_column}      = 7;
        $option{y_last_column} = 7;
        $option{y_title}       = "Di/Dn";
        $option{y_scale_unit}  = 0.1;
        $option{Top} += $option{Height} + 12.75;

        $excel_obj->draw_y( $sheet_name, \%option );
    }
}

#----------------------------------------------------------#
# POST Processing
#----------------------------------------------------------#
# add time stamp to "summary" sheet
$excel_obj->time_stamp("summary") if $time_stamp;

# add an index sheet
$excel_obj->add_index_sheet if $add_index_sheet;

print "$outfile has been generated.\n";

$stopwatch->end_message();
exit;

__END__

    'help|?'   => \$help,
    'man'      => \$man,
    'infile=s' => \$infile,
    'outfile=s' => \$outfile,
    
=head1 NAME

    three_chart_factory.pl - Use Win32::OLE to automate Excel chart

=head1 SYNOPSIS

    three_chart_factory.pl [options]
     Options:
       --help            brief help message
       --man             full documentation
       --infile          input file name (full path)
       --outfile         output file name
       

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do someting
useful with the contents thereof.

=cut
