#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB::Excel;
use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $infile          = '';
my $outfile         = '';
my $jc_correction   = 0;
my $time_stamp      = 1;
my $add_index_sheet = 1;

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

my $excel_obj = AlignDB::Excel->new( infile => $infile, );
if ($outfile) {
    $excel_obj->outfile($outfile);
}

#----------------------------------------------------------#
# START
#----------------------------------------------------------#
# jc
$excel_obj->jc_correction if $jc_correction;

#----------------------------------------------------------#
# draw charts section
#----------------------------------------------------------#
my @sheet_names = @{ $excel_obj->sheet_names };
{

    #----------------------------#
    # worksheet -- distance_
    #----------------------------#
    my @sheets = grep {/^distance/} @sheet_names;
    foreach (@sheets) {
        my $sheet_name = $_;
        my %option     = (
            chart_serial  => 1,
            x_column      => 1,
            y_column      => 3,
            y_last_column => 6,
            first_row     => 3,
            last_row      => 8,
            x_max_scale   => 5,
            x_title       => "Distance to indels (d1)",
            y_title       => "Nucleotide diversity",
            Height        => 200,
            Width         => 320,
            Top           => 12.75,
            Left          => 520,
        );
        $excel_obj->draw_y( $sheet_name, \%option );

        # chart 2
        $option{chart_serial}++;
        $option{y_column}      = 8;
        $option{y_last_column} = 8;
        $option{y_title}       = "Di/Dn";
        $option{Top} += $option{Height} + 12.75;
        $excel_obj->draw_y( $sheet_name, \%option );
    }
}

{

    #----------------------------#
    # worksheet -- combined_pigccv
    #----------------------------#
    my $sheet_name = 'combined_pigccv';
    my %option     = (
        chart_serial => 1,
        x_column     => 1,
        y_column     => 2,
        first_row    => 3,
        last_row     => 33,
        x_max_scale  => 30,
        x_title      => "Distance to indels (d1)",
        y_title      => "Nucleotide diversity",
        Height       => 200,
        Width        => 320,
        Top          => 12.75,
        Left         => 520,
    );
    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 2
    $option{chart_serial}++;
    $option{y_column} = 4;
    $option{y_title}  = "GC proportion";
    $option{Top} += $option{Height} + 12.75;
    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 3
    $option{chart_serial}++;
    $option{y_column} = 6;
    $option{y_title}  = "CV";
    $option{Top} += $option{Height} + 12.75;
    $excel_obj->draw_y( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- combined_pure_coding
    #----------------------------#
    $sheet_name           = 'combined_pure_coding';
    $option{chart_serial} = 1;
    $option{y_column}     = 2;
    $option{y_title}      = "Nucleotide diversity";
    $option{Top}          = 12.75;
    $excel_obj->draw_y( $sheet_name, \%option );
}

{

    #----------------------------#
    # worksheet -- pigccv_freq_
    #----------------------------#
    my @sheets = grep {/^pigccv_freq/} @sheet_names;
    foreach (@sheets) {
        my $sheet_name = $_;
        my %option     = (
            chart_serial => 1,
            x_column     => 1,
            y_column     => 2,
            first_row    => 3,
            last_row     => 8,
            x_max_scale  => 5,
            x_title      => "Distance to indels (d1)",
            y_title      => "Nucleotide diversity",
            Height       => 200,
            Width        => 320,
            Top          => 12.75,
            Left         => 520,
        );
        $excel_obj->draw_y( $sheet_name, \%option );

        # chart 2
        $option{chart_serial}++;
        $option{y_column} = 4;
        $option{y_title}  = "GC proportion";
        $option{Top} += $option{Height} + 12.75;
        $excel_obj->draw_y( $sheet_name, \%option );

        # chart 3
        $option{chart_serial}++;
        $option{y_column} = 6;
        $option{y_title}  = "CV";
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
    
=head1 NAME

    multi_chart_factory.pl - Use Win32::OLE to automate Excel chart

=head1 SYNOPSIS

    multi_chart_factory.pl [options]
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
