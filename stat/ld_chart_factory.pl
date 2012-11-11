#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB::Excel;
use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

my $infile        = '';
my $outfile       = '';
my $jc_correction = 0;

my %replace;

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'    => \$help,
    'man'       => \$man,
    'infile=s'  => \$infile,
    'outfile=s' => \$outfile,
    'jc=s'      => \$jc_correction,
    'replace=s' => \%replace,
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
my @sheet_names = @{ $excel_obj->sheet_names };

{

    #----------------------------#
    # worksheet -- ld
    #----------------------------#
    my @sheets = grep {/^ld/} @sheet_names;
    foreach (@sheets) {
        my $sheet_name = $_;
        my %option     = (
            chart_serial => 1,
            x_column     => 1,
            y_column     => 2,
            y2_column    => 4,
            first_row    => 3,
            last_row     => 13,
            x_max_scale  => 10,
            x_title      => "Distance to indels (d1)",
            y_title      => "r",
            y2_title     => "Dprime",
            Height       => 200,
            Width        => 260,
            Top          => 14.25,
            Left         => 650,
        );
        $excel_obj->draw_2y( $sheet_name, \%option );

        # chart 2
        $option{chart_serial}++;
        $option{y_column}  = 3;
        $option{y_title}   = "r**2";
        $option{y2_column} = 5;
        $option{y2_title}  = "|Dprime|";
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_2y( $sheet_name, \%option );
    }
}

{

    #----------------------------#
    # worksheet -- sub_ld
    #----------------------------#
    my @sheets = grep {/^sub_ld/} @sheet_names;
    foreach (@sheets) {
        my $sheet_name = $_;
        my %option     = (
            chart_serial => 1,
            x_column     => 1,
            y_column     => 2,
            y2_column    => 4,
            first_row    => 3,
            last_row     => 13,
            x_max_scale  => 10,
            x_title      => "Distance to indels (d1)",
            y_title      => "indel-group r**2",
            y2_title     => "indel-group |Dprime|",
            Height       => 200,
            Width        => 260,
            Top          => 14.25,
            Left         => 650,
        );
        $excel_obj->draw_2y( $sheet_name, \%option );

        # chart 2
        $option{chart_serial}++;
        $option{y_column}  = 3;
        $option{y_title}   = "nonindel-group r**2",
        $option{y2_column} = 5;
        $option{y2_title}  = "nonindel-group |Dprime|",
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_2y( $sheet_name, \%option );
    }
}

{

    #----------------------------#
    # worksheet -- snp_ld
    #----------------------------#
    my @sheets = grep {/^snp_ld/} @sheet_names;
    foreach (@sheets) {
        my $sheet_name = $_;
        my %option     = (
            chart_serial => 1,
            x_column     => 1,
            y_column     => 2,
            y2_column    => 4,
            first_row    => 3,
            last_row     => 13,
            x_max_scale  => 10,
            x_title      => "Distance to indels (d1)",
            y_title      => "nearest indel r**2",
            y2_title     => "nearest indel |Dprime|",
            Height       => 200,
            Width        => 260,
            Top          => 14.25,
            Left         => 650,
        );
        $excel_obj->draw_2y( $sheet_name, \%option );

        # chart 2
        $option{chart_serial}++;
        $option{y_column}  = 3;
        $option{y_title}   = "nearest snp r**2",
        $option{y2_column} = 5;
        $option{y2_title}  = "nearest snp |Dprime|",
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_2y( $sheet_name, \%option );
    }
}

{

    #----------------------------#
    # worksheet -- snps_ld
    #----------------------------#
    my @sheets = grep {/^snps_ld/} @sheet_names;
    foreach (@sheets) {
        my $sheet_name = $_;
        my %option     = (
            chart_serial => 1,
            x_column     => 1,
            y_column     => 2,
            y2_column    => 4,
            first_row    => 3,
            last_row     => 13,
            x_max_scale  => 10,
            x_title      => "Distance to indels (d1)",
            y_title      => "nearest indel r**2",
            y2_title     => "nearest indel |Dprime|",
            Height       => 200,
            Width        => 260,
            Top          => 14.25,
            Left         => 650,
        );
        $excel_obj->draw_2y( $sheet_name, \%option );

        # chart 2
        $option{chart_serial}++;
        $option{y_column}  = 3;
        $option{y_title}   = "near snps r**2",
        $option{y2_column} = 5;
        $option{y2_title}  = "near snps |Dprime|",
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_2y( $sheet_name, \%option );
    }
}

#----------------------------------------------------------#
# POST Processing
#----------------------------------------------------------#
# add time stamp to "basic" sheet
$excel_obj->time_stamp("basic");

# add an index sheet
$excel_obj->add_index_sheet;

print "$outfile has been generated.\n";

$stopwatch->end_message;
exit;

__END__
    
=head1 NAME

    ld_chart_factory.pl - Use Win32::OLE to automate Excel chart

=head1 SYNOPSIS

    ld_chart_factory.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --infile            input file name (full path)
        --outfile           output file name
        --jc                Jukes & Cantor correction
        --replace           replace text when charting
                            --replace diversity=divergence

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
