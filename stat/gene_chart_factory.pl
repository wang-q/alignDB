#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

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

my $infile  = '';
my $outfile = '';

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
    # worksheet -- gene_ess, exon_ess, exon_gc
    #----------------------------#
    my @sheets = qw{
        coding_all
        gene_ess_ess gene_ess_non_ess
        gene_pure_ess gene_pure_non_ess
        gene_5_ess gene_5_non_ess
        gene_3_ess gene_3_non_ess
        exon_ess_ess exon_ess_non_ess
        exon_gc_1 exon_gc_2
        exon_gc_3 exon_gc_4
        exon_gc_all
        gene_all_pure coding_all_pure
        gene_D_1 gene_D_2
        gene_D_3 gene_D_4
        gene_D_all gene_D_null
        exon_D_1 exon_D_2
        exon_D_3 exon_D_4
        exon_D_all exon_D_null
        coding_quan_1 coding_quan_2
        coding_quan_3 coding_quan_4
    };

    foreach (@sheets) {
        my $sheet_name = $_;
        my $x_title
            = /gene/ ? "Distance to gene"
            : /exon/ ? "Distance to exon"
            :          "Distance to coding";
        my %option = (
            chart_serial => 1,
            x_column     => 1,
            y_column     => 2,
            first_row    => 2,
            last_row     => 17,
            x_min_scale  => -5,
            x_max_scale  => 10,
            x_title      => $x_title,
            y_title      => "Nucleotide diversity",
            cross        => -5,
            Height       => 200,
            Width        => 320,
            Top          => 12.75,
            Left         => 450,
        );

        $excel_obj->draw_y( $sheet_name, \%option );

        # chart 2
        $option{chart_serial}++;
        $option{y_column} = 3;
        $option{y_title}  = "Indel per 100 bp";
        $option{Top} += $option{Height} + 12.75;
        $excel_obj->draw_y( $sheet_name, \%option );

        # chart 3
        $option{chart_serial}++;
        $option{y_column} = 4;
        $option{y_title}  = "GC proportion";
        $option{Top} += $option{Height} + 12.75;
        $excel_obj->draw_y( $sheet_name, \%option );

        # chart 4
        $option{chart_serial}++;
        $option{y_column} = 5;
        $option{y_title}  = "Window CV";
        $option{Top} += $option{Height} + 12.75;
        $excel_obj->draw_y( $sheet_name, \%option );

        # chart 4
        $option{chart_serial}++;
        $option{y_column} = 6;
        $option{y_title}  = "Repeats proportion";
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

$stopwatch->end_message;
exit;

__END__
    
=head1 NAME

    gene_chart_factory.pl - Use Win32::OLE to automate Excel chart

=head1 SYNOPSIS

    gene_chart_factory.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --infile            input file name (full path)
        --outfile           output file name

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
