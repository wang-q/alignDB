#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB::Excel;
use AlignDB::Stopwatch;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $Config = Config::Tiny->new;
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# stat parameter
my $jc_correction = $Config->{stat}{jc_correction};

my $infile  = '';
my $outfile = '';

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

if (! -e $infile) {
    warn "[$infile] doesn't exist.\n";
    exit;
}

#----------------------------------------------------------#
# Init section
#----------------------------------------------------------#
if (! -e $infile) {
    warn "[$infile] doesn't exist.\n";
    exit;
}

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
    my %option = (
        chart_serial => 1,
        x_column     => 1,
        y_column     => 2,
        first_row    => 3,
        last_row     => 18,
        x_max_scale  => 15,
        x_title      => "Distance to indels (D1)",
        y_title      => "Nucleotide diversity",
        Height       => 200,
        Width        => 260,
        Top          => 14.25,
        Left         => 520,
    );

    #----------------------------#
    # worksheet -- align_coding
    #----------------------------#
    my @coding_levels = ( 1 .. 9 );

    foreach (@coding_levels) {
        my $sheet_name = "align_coding_$_";

        $excel_obj->draw_y( $sheet_name, \%option );
    }

    #----------------------------#
    # worksheet -- align_repeat
    #----------------------------#
    my @repeat_levels = ( 1 .. 9 );

    foreach (@repeat_levels) {
        my $sheet_name = "align_repeat_$_";

        $excel_obj->draw_y( $sheet_name, \%option );
    }
}

{

    #----------------------------#
    # worksheet -- density_snp
    #----------------------------#
    my $sheet_name = 'density_snp';
    my %option     = (
        chart_serial  => 1,
        x_column      => 1,
        y_column      => 2,
        y_last_column => 7,
        first_row     => 3,
        last_row      => 33,
        x_max_scale   => 30,
        x_title       => "Indel density (d2)",
        y_title       => "Proportion of substitutions",
        Height        => 300,
        Width         => 390,
        Top           => 14.25,
        Left          => 520,
    );

    $excel_obj->draw_y( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- distance_snp
    #----------------------------#
    $sheet_name          = 'distance_snp';
    $option{last_row}    = 18;
    $option{x_max_scale} = 15;
    $option{x_title}     = "Distance to indels (d1)";

    $excel_obj->draw_y( $sheet_name, \%option );
}

#----------------------------------------------------------#
# POST Processing
#----------------------------------------------------------#
# add time stamp to "summary" sheet
$excel_obj->time_stamp("basic");

# add an index sheet
$excel_obj->add_index_sheet;

print "$outfile has been generated.\n";

$stopwatch->end_message;
exit;

__END__
    
=head1 NAME

    common_chart_factory.pl - Use Win32::OLE to automate Excel chart

=head1 SYNOPSIS

    common_chart_factory.pl [options]
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
