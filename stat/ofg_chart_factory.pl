#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Config::Tiny;
use YAML::Syck qw(Dump Load DumpFile LoadFile);

use Set::Scalar;
use List::MoreUtils qw( first_index);

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
my @sheet_names = @{ $excel_obj->sheet_names };
{

    #----------------------------#
    # worksheet -- gene_ess, exon_ess, exon_gc
    #----------------------------#
    my @sheets = grep {/^ofg/} @sheet_names;

    foreach (@sheets) {
        my $sheet_name = $_;

        my $x_column = 1;
        my $values   = $excel_obj->get_column( $sheet_name, $x_column );
        my $set      = Set::Scalar->new( @{$values} );

        my %option = (
            chart_serial => 1,
            x_column     => $x_column,
            Height       => 200,
            Width        => 260,
            Top          => 12.75,
            Left         => 800,
        );

        if ( $set->has(-10) ) {
            my $idx = first_index { $_ == -10 } @{$values};
            $option{first_row}   = $idx + 2;
            $option{last_row}    = $idx + 32;
            $option{x_min_scale} = -10;
            $option{x_max_scale} = 20;
            $option{cross}       = 0;
        }
        elsif ( $set->has(0) ) {
            my $idx = first_index { $_ == 0 } @{$values};
            $option{first_row}   = $idx + 2;
            $option{last_row}    = $idx + 32;
            $option{x_min_scale} = 0;
            $option{x_max_scale} = 20;
            $option{cross}       = 0;
        }
        elsif ( $set->has(1) ) {
            my $idx = first_index { $_ == 1 } @{$values};
            $option{first_row}   = $idx + 2;
            $option{last_row}    = $idx + 32;
            $option{x_min_scale} = 0;
            $option{x_max_scale} = 10;
        }
        else {
            print "X column errors\n";
            next;
        }

        $option{y_column} = 2;
        $option{x_title}  = "Distance to ofg";
        $option{y_title}  = "Nucleotide diversity";
        $excel_obj->draw_y( $sheet_name, \%option );

        # chart 2
        $option{chart_serial}++;
        $option{y_column} = 4;
        $option{y_title}  = "Indel per 100 bp";
        $option{Top} += $option{Height} + 12.75;
        $excel_obj->draw_y( $sheet_name, \%option );

        # chart 3
        $option{chart_serial}++;
        $option{y_column} = 6;
        $option{y_title}  = "GC proportion";
        $option{Top} += $option{Height} + 12.75;
        $excel_obj->draw_y( $sheet_name, \%option );

        # chart 4
        $option{chart_serial}++;
        $option{y_column} = 8;
        $option{y_title}  = "CV";
        $option{Top} += $option{Height} + 12.75;
        $excel_obj->draw_y( $sheet_name, \%option );

        # chart 5
        $option{chart_serial}++;
        $option{y_column} = 10;
        $option{y_title}  = "Repeats proportion";
        $option{Top} += $option{Height} + 12.75;
        $excel_obj->draw_y( $sheet_name, \%option );

        if ( $set->has(-90) ) {
            my $idx = first_index { $_ == -90 } @{$values};
            $option{first_row}   = $idx + 2;
            $option{last_row}    = $idx + 22;
            $option{x_min_scale} = -90;
            $option{x_max_scale} = -70;
            $option{cross}       = -90;

            # chart 6
            $option{chart_serial}++;
            $option{y_column} = 2;
            $option{x_title}  = "Distance to ofg";
            $option{y_title}  = "Nucleotide diversity";
            $option{Top}      = 12.75;
            $option{Left}     = 1100;
            $excel_obj->draw_y( $sheet_name, \%option );

            # chart 7
            $option{chart_serial}++;
            $option{y_column} = 4;
            $option{y_title}  = "Indel per 100 bp";
            $option{Top} += $option{Height} + 12.75;
            $excel_obj->draw_y( $sheet_name, \%option );

            # chart 8
            $option{chart_serial}++;
            $option{y_column} = 6;
            $option{y_title}  = "GC proportion";
            $option{Top} += $option{Height} + 12.75;
            $excel_obj->draw_y( $sheet_name, \%option );

            # chart 9
            $option{chart_serial}++;
            $option{y_column} = 8;
            $option{y_title}  = "CV";
            $option{Top} += $option{Height} + 12.75;
            $excel_obj->draw_y( $sheet_name, \%option );

            # chart 10
            $option{chart_serial}++;
            $option{y_column} = 10;
            $option{y_title}  = "Repeats proportion";
            $option{Top} += $option{Height} + 12.75;
            $excel_obj->draw_y( $sheet_name, \%option );
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
