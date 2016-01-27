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

if ( !-e $infile ) {
    warn "[$infile] doesn't exist.\n";
    exit;
}

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
