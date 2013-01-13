#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML qw(Dump Load DumpFile LoadFile);

use AlignDB::Excel;
use AlignDB::Stopwatch;

use FindBin;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# stat parameter
my $jc_correction   = $Config->{stat}{jc_correction};
my $time_stamp      = $Config->{stat}{time_stamp};
my $add_index_sheet = $Config->{stat}{add_index_sheet};

my $add_trend = 0;

my $infile  = '';
my $outfile = '';

my %replace;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'            => \$help,
    'man'               => \$man,
    'infile=s'          => \$infile,
    'outfile=s'         => \$outfile,
    'jc=s'              => \$jc_correction,
    'time_stamp=s'      => \$time_stamp,
    'add_index_sheet=s' => \$add_index_sheet,
    'replace=s'         => \%replace,
    'add_trend=s'       => \$add_trend,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Init section
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;
$stopwatch->start_message("Processing $infile...");

my $excel_obj;
if ($outfile) {
    $excel_obj = AlignDB::Excel->new(
        infile  => $infile,
        outfile => $outfile,
        replace => \%replace,
    );
}
else {
    $excel_obj = AlignDB::Excel->new(
        infile  => $infile,
        replace => \%replace,
    );
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
    # worksheet -- distance_to_trough
    #----------------------------#
    my $sheet_name = 'distance_to_trough';
    my %option     = (
        chart_serial => 1,
        x_column     => 1,
        y_column     => 2,
        first_row    => 2,
        last_row     => 17,
        x_max_scale  => 15,
        x_title      => "Distance to GC trough",
        y_title      => "Nucleotide diversity",
        Height       => 200,
        Width        => 260,
        Top          => 14.25,
        Left         => 550,
    );
    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 2
    $option{chart_serial}++;
    $option{y_column} = 4;
    $option{y_title}  = "Indel per 100 bp";
    $option{Top} += $option{Height} + 14.25;
    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 3
    $option{chart_serial}++;
    $option{y_column} = 6;
    $option{y_title}  = "Window CV";
    $option{Top} += $option{Height} + 14.25;
    $excel_obj->draw_y( $sheet_name, \%option );
}

{
    #----------------------------#
    # worksheet -- wave_length
    #----------------------------#
    my $sheet_name = 'wave_length';
    my %option     = (
        chart_serial => 1,
        x_column     => 1,
        y_column     => 2,
        first_row    => 2,
        last_row     => 32,
        x_title      => "Wave length",
        y_title      => "Nucleotide diversity",
        Height       => 200,
        Width        => 260,
        Top          => 14.25,
        Left         => 550,
    );
    $option{x_min_scale}  = 5;
    $option{x_max_scale}  = 35;
    $option{x_scale_unit} = 5;
    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 2
    $option{chart_serial}++;
    $option{y_column} = 4;
    $option{y_title}  = "Indel per 100 bp";
    $option{Top} += $option{Height} + 14.25;
    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 3
    $option{chart_serial}++;
    $option{y_column} = 6;
    $option{y_title}  = "Window CV";
    $option{Top} += $option{Height} + 14.25;
    $excel_obj->draw_y( $sheet_name, \%option );
}

{
    #----------------------------#
    # worksheet -- amplitude
    #----------------------------#
    my $sheet_name = 'amplitude';
    my %option     = (
        chart_serial => 1,
        x_column     => 1,
        y_column     => 2,
        first_row    => 2,
        last_row     => 32,
        x_title      => "Amplitude",
        y_title      => "Nucleotide diversity",
        Height       => 200,
        Width        => 260,
        Top          => 14.25,
        Left         => 550,
    );
    $option{x_min_scale}  = 10;
    $option{x_max_scale}  = 40;
    $option{x_scale_unit} = 5;
    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 2
    $option{chart_serial}++;
    $option{y_column} = 4;
    $option{y_title}  = "Indel per 100 bp";
    $option{Top} += $option{Height} + 14.25;
    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 3
    $option{chart_serial}++;
    $option{y_column} = 6;
    $option{y_title}  = "Window CV";
    $option{Top} += $option{Height} + 14.25;
    $excel_obj->draw_y( $sheet_name, \%option );
}

{
    #----------------------------#
    # worksheet -- trough_gc
    #----------------------------#
    my $sheet_name = 'trough_gc';
    my %option     = (
        chart_serial => 1,
        x_column     => 1,
        y_column     => 2,
        first_row    => 2,
        last_row     => 32,
        x_title      => "Trough GC",
        y_title      => "Nucleotide diversity",
        Height       => 200,
        Width        => 260,
        Top          => 14.25,
        Left         => 550,
    );
    $option{x_title}      = "Trough GC";
    $option{y_title}      = "Nucleotide diversity";
    $option{x_scale_unit} = 5;
    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 2
    $option{chart_serial}++;
    $option{y_column} = 4;
    $option{y_title}  = "Indel per 100 bp";
    $option{Top} += $option{Height} + 14.25;
    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 3
    $option{chart_serial}++;
    $option{y_column} = 6;
    $option{y_title}  = "Window CV";
    $option{Top} += $option{Height} + 14.25;
    $excel_obj->draw_y( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- gradient
    #----------------------------#
    $sheet_name           = 'gradient';
    $option{chart_serial} = 1;
    $option{y_column}     = 2;
    $option{x_title}      = "Gradient";
    $option{y_title}      = "Nucleotide diversity";
    $option{Top}          = 14.25;
    $option{x_max_scale}  = 15;
    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 2
    $option{chart_serial}++;
    $option{y_column} = 4;
    $option{y_title}  = "Indel per 100 bp";
    $option{Top} += $option{Height} + 14.25;
    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 3
    $option{chart_serial}++;
    $option{y_column} = 6;
    $option{y_title}  = "Window CV";
    $option{Top} += $option{Height} + 14.25;
    $excel_obj->draw_y( $sheet_name, \%option );
}

{

    #----------------------------#
    # worksheet -- d_wave_length_series
    #----------------------------#
    my $sheet_name = 'd_wave_length_series';
    my %option     = (
        chart_serial   => 1,
        x_title        => "Distance to GC trough",
        y_title        => "Indel per 100 bp",
        Height         => 300,
        Width          => 390,
        Top            => 14.25 * 17,
        Left           => 360,
        section_top    => 2,
        section_end    => 13,
        section_length => 12,
        x_orientation  => 0,
    );
    $option{group_name} = [1 .. 4] ;
    $excel_obj->draw_dd( $sheet_name, \%option );
    
    #----------------------------#
    # worksheet -- d_amplitude_series
    #----------------------------#
    $sheet_name = 'd_amplitude_series';
    $option{group_name} = [1 .. 4] ;
    $excel_obj->draw_dd( $sheet_name, \%option );
    
    #----------------------------#
    # worksheet -- d_gradient_series
    #----------------------------#
    $sheet_name = 'd_gradient_series';
    $option{group_name} = [1 .. 4] ;
    $excel_obj->draw_dd( $sheet_name, \%option );
    
    #----------------------------#
    # worksheet -- dtrough_gc_series
    #----------------------------#
    $sheet_name = 'd_trough_gc_series';
    $option{group_name} = [1 .. 4] ;
    $excel_obj->draw_dd( $sheet_name, \%option );
    
    #----------------------------#
    # worksheet -- dgc_series
    #----------------------------#
    $sheet_name = 'd_gc_series';
    $option{group_name} = [1 .. 6] ;
    $excel_obj->draw_dd( $sheet_name, \%option );
}

{

    #----------------------------#
    # worksheet -- extreme_amplitude_group
    #----------------------------#
    # select worksheet by name
    my $sheet_name = 'extreme_amplitude_group';
    my @group_name = (
        "Descend 0-0.15 & Ascend 0-0.15",
        "Descend 0.15-1 & Ascend 0-0.15",
        "Descend 0-0.15 & Ascend 0.15-1",
        "Descend 0.15-1 & Ascend 0.15-1",
    );
    my %option = (
        chart_serial   => 1,
        x_title        => "Distance to GC trough",
        y_title        => "Indel per 100 bp",
        Height         => 300,
        Width          => 390,
        Top            => 14.25 * 17,
        Left           => 360,
        group_name     => \@group_name,
        section_top    => 2,
        section_end    => 15,
        section_length => 14,
        x_orientation  => 0,
    );

    $excel_obj->draw_dd( $sheet_name, \%option );
}


{

    #----------------------------#
    # worksheet -- segment_gc_indel_x
    #----------------------------#
    foreach ( 'A', 0 .. 9, 10, 20, 30, 40, 50 ) {
        my $sheet_name = 'segment_gc_indel_' . $_;
        my %option     = (
            chart_serial => 1,
            x_column     => 2,
            y_column     => 3,
            x_title      => "GC proportion",
            y_title      => "Nucleotide diversity",
            Height       => 200,
            Width        => 260,
            Top          => 14.25,
            Left         => 650,
            without_line => 1,
            marker_size  => 5,
            add_trend    => $add_trend,
        );
        $excel_obj->draw_xy( $sheet_name, \%option );
        $excel_obj->linear_fit( $sheet_name, \%option ) if $_ eq 3;

        # chart 2
        $option{chart_serial}++;
        $option{y_column} = 4;
        $option{y_title}  = "Indel per 100bp";
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_xy( $sheet_name, \%option );
        $excel_obj->linear_fit( $sheet_name, \%option ) if $_ eq 3;

        # chart 3
        $option{chart_serial}++;
        $option{y_column} = 5;
        $option{y_title}  = "Segment CV";
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_xy( $sheet_name, \%option );
        $excel_obj->linear_fit( $sheet_name, \%option ) if $_ eq 3;

        # chart 4
        $option{chart_serial}++;
        $option{y_column} = 6;
        $option{y_title}  = "Coding proportion";
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_xy( $sheet_name, \%option );
    }

    #----------------------------#
    # worksheet -- segment_std_indel_x
    #----------------------------#
    foreach ( 'A', 0 .. 9, 10, 20, 30, 40, 50 ) {
        my $sheet_name = 'segment_std_indel_' . $_;
        my %option     = (
            chart_serial => 1,
            x_column     => 2,
            y_column     => 3,
            x_title      => "Segment std",
            y_title      => "Nucleotide diversity",
            Height       => 200,
            Width        => 260,
            Top          => 14.25,
            Left         => 650,
            without_line => 1,
            marker_size  => 5,
            add_trend    => $add_trend,
        );
        $excel_obj->draw_xy( $sheet_name, \%option );
        $excel_obj->linear_fit( $sheet_name, \%option ) if $_ eq 3;

        # chart 2
        $option{chart_serial}++;
        $option{y_column} = 4;
        $option{y_title}  = "Indel per 100bp";
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_xy( $sheet_name, \%option );
        $excel_obj->linear_fit( $sheet_name, \%option ) if $_ eq 3;

        # chart 3
        $option{chart_serial}++;
        $option{y_column} = 5;
        $option{y_title}  = "GC proportion";
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_xy( $sheet_name, \%option );
        $excel_obj->linear_fit( $sheet_name, \%option ) if $_ eq 3;

        # chart 4
        $option{chart_serial}++;
        $option{y_column} = 6;
        $option{y_title}  = "Coding proportion";
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_xy( $sheet_name, \%option );
    }

    #----------------------------#
    # worksheet -- segment_cv_indel_x
    #----------------------------#
    foreach ( 'A', 0 .. 9, 10, 20, 30, 40, 50 ) {
        my $sheet_name = 'segment_cv_indel_' . $_;
        my %option     = (
            chart_serial => 1,
            x_column     => 2,
            y_column     => 3,
            x_title      => "Segment CV",
            y_title      => "Nucleotide diversity",
            Height       => 200,
            Width        => 260,
            Top          => 14.25,
            Left         => 650,
            without_line => 1,
            marker_size  => 5,
            add_trend    => $add_trend,
        );
        $excel_obj->draw_xy( $sheet_name, \%option );
        $excel_obj->linear_fit( $sheet_name, \%option ) if $_ eq 3;

        # chart 2
        $option{chart_serial}++;
        $option{y_column} = 4;
        $option{y_title}  = "Indel per 100bp";
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_xy( $sheet_name, \%option );
        $excel_obj->linear_fit( $sheet_name, \%option ) if $_ eq 3;

        # chart 3
        $option{chart_serial}++;
        $option{y_column} = 5;
        $option{y_title}  = "GC proportion";
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_xy( $sheet_name, \%option );
        $excel_obj->linear_fit( $sheet_name, \%option ) if $_ eq 3;

        # chart 4
        $option{chart_serial}++;
        $option{y_column} = 6;
        $option{y_title}  = "Coding proportion";
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_xy( $sheet_name, \%option );

        # chart 5
        $option{chart_serial}++;
        $option{y_column} = 10;
        $option{y_title}  = "GC Range";
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_xy( $sheet_name, \%option );
    }

    #----------------------------#
    # worksheet -- segment_mdcw_indel_x
    #----------------------------#
    foreach ( 'A', 0 .. 9, 10, 20, 30, 40, 50 ) {
        my $sheet_name = 'segment_mdcw_indel_' . $_;
        my %option     = (
            chart_serial => 1,
            x_column     => 2,
            y_column     => 3,
            x_title      => "Segment mdcw",
            y_title      => "Nucleotide diversity",
            Height       => 200,
            Width        => 260,
            Top          => 14.25,
            Left         => 650,
            without_line => 1,
            marker_size  => 5,
            add_trend    => $add_trend,
        );
        $excel_obj->draw_xy( $sheet_name, \%option );
        $excel_obj->linear_fit( $sheet_name, \%option ) if $_ eq 3;

        # chart 2
        $option{chart_serial}++;
        $option{y_column} = 4;
        $option{y_title}  = "Indel per 100bp";
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_xy( $sheet_name, \%option );
        $excel_obj->linear_fit( $sheet_name, \%option ) if $_ eq 3;

        # chart 3
        $option{chart_serial}++;
        $option{y_column} = 5;
        $option{y_title}  = "GC proportion";
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_xy( $sheet_name, \%option );
        $excel_obj->linear_fit( $sheet_name, \%option ) if $_ eq 3;

        # chart 4
        $option{chart_serial}++;
        $option{y_column} = 6;
        $option{y_title}  = "Coding proportion";
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_xy( $sheet_name, \%option );
    }

    #----------------------------#
    # worksheet -- segment_coding_indel_x
    #----------------------------#
    foreach ( 'A', 0 .. 9, 10, 20, 30, 40, 50 ) {
        my $sheet_name = 'segment_coding_indel_' . $_;
        my %option     = (
            chart_serial => 1,
            x_column     => 2,
            y_column     => 3,
            x_title      => "Segment coding",
            y_title      => "Nucleotide diversity",
            Height       => 200,
            Width        => 260,
            Top          => 14.25,
            Left         => 600,
            without_line => 1,
            marker_size  => 5,
        );
        $excel_obj->draw_xy( $sheet_name, \%option );

        # chart 2
        $option{chart_serial}++;
        $option{y_column} = 4;
        $option{y_title}  = "Indel per 100 bp";
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_xy( $sheet_name, \%option );

        # chart 3
        $option{chart_serial}++;
        $option{y_column} = 5;
        $option{y_title}  = "GC proportion";
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_xy( $sheet_name, \%option );

        # chart 4
        $option{chart_serial}++;
        $option{y_column} = 6;
        $option{y_title}  = "Segment CV";
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_xy( $sheet_name, \%option );
    }

    #----------------------------#
    # worksheet -- segment_gc_indel_cr_3_xx
    #----------------------------#
    for my $i (3) {
        for my $j (qw{ 10 01 00 }) {
            my $sheet_name = 'segment_gc_indel_cr_' . $i . '_' . $j;
            my %option     = (
                chart_serial => 1,
                x_column     => 2,
                y_column     => 3,
                x_title      => "GC proportion",
                y_title      => "Nucleotide diversity",
                Height       => 200,
                Width        => 260,
                Top          => 14.25,
                Left         => 600,
                without_line => 1,
                marker_size  => 5,
            );
            $excel_obj->draw_xy( $sheet_name, \%option );

            # chart 2
            $option{chart_serial}++;
            $option{y_column} = 4;
            $option{y_title}  = "Indel per 100 bp";
            $option{Top} += $option{Height} + 14.25;
            $excel_obj->draw_xy( $sheet_name, \%option );

            # chart 3
            $option{chart_serial}++;
            $option{y_column} = 5;
            $option{y_title}  = "Segment CV";
            $option{Top} += $option{Height} + 14.25;
            $excel_obj->draw_xy( $sheet_name, \%option );

            # chart 4
            $option{chart_serial}++;
            $option{y_column} = 6;
            $option{y_title}  = "Coding proportion";
            $option{Top} += $option{Height} + 14.25;
            $excel_obj->draw_xy( $sheet_name, \%option );
        }
    }

    #----------------------------#
    # worksheet -- segment_cv_indel_cr_3_xx
    #----------------------------#
    for my $i (3) {
        for my $j (qw{ 10 01 00 }) {
            my $sheet_name = 'segment_cv_indel_cr_' . $i . '_' . $j;
            my %option     = (
                chart_serial => 1,
                x_column     => 2,
                y_column     => 3,
                x_title      => "Segment CV",
                y_title      => "Nucleotide diversity",
                Height       => 200,
                Width        => 260,
                Top          => 14.25,
                Left         => 600,
                without_line => 1,
                marker_size  => 5,
            );
            $excel_obj->draw_xy( $sheet_name, \%option );

            # chart 2
            $option{chart_serial}++;
            $option{y_column} = 4;
            $option{y_title}  = "Indel per 100 bp";
            $option{Top} += $option{Height} + 14.25;
            $excel_obj->draw_xy( $sheet_name, \%option );

            # chart 3
            $option{chart_serial}++;
            $option{y_column} = 5;
            $option{y_title}  = "GC proportion";
            $option{Top} += $option{Height} + 14.25;
            $excel_obj->draw_xy( $sheet_name, \%option );

            # chart 4
            $option{chart_serial}++;
            $option{y_column} = 6;
            $option{y_title}  = "Coding proportion";
            $option{Top} += $option{Height} + 14.25;
            $excel_obj->draw_xy( $sheet_name, \%option );
        }
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

$stopwatch->end_message;
exit;

__END__
    
=head1 NAME

    gc_chart_factory.pl - Use Win32::OLE to automate Excel chart

=head1 SYNOPSIS

    gc_chart_factory.pl [options]
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
