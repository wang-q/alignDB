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
my $jc_correction   = $Config->{stat}{jc_correction};
my $time_stamp      = $Config->{stat}{time_stamp};
my $add_index_sheet = $Config->{stat}{add_index_sheet};

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
    # worksheet -- combined_distance
    #----------------------------#
    my @sheets = qw{
        combined_distance distance_coding distance_non_coding
        combined_density density_coding density_non_coding
    };
    for my $sheet_name (@sheets) {
        my %option = (
            chart_serial => 1,
            x_column     => 1,
            y_column     => 2,
            first_row    => 3,
            last_row     => 23,
            x_max_scale  => 20,
            x_title      => "Distance to indels (d1)",
            y_title      => "Nucleotide diversity",
            Height       => 200,
            Width        => 260,
            Top          => 14.25,
            Left         => 520,
        );
        if ( $sheet_name =~ /density/ ) {
            $option{x_title}     = "Reciprocal of indel density (d2)";
            $option{last_row}    = 43;
            $option{x_max_scale} = 40;
        }
        $excel_obj->draw_y( $sheet_name, \%option );

        # chart 2
        $option{chart_serial}++;
        $option{y_column} = 4;
        $option{y_title}  = "GC proportion";
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_y( $sheet_name, \%option );

        # chart 3
        $option{chart_serial}++;
        $option{y_column} = 6;
        $option{y_title}  = "CV";
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_y( $sheet_name, \%option );

        # chart 4
        $option{chart_serial}++;
        $option{y_column} = 4;
        $option{y_title}  = "GC proportion";
        $option{y2_column} = 6;
        $option{y2_title}  = "CV";
        $option{Top} += $option{Height} + 14.25;
        $excel_obj->draw_2y( $sheet_name, \%option );
    }
}

{

    #----------------------------#
    # worksheet -- combined_distance
    #----------------------------#
    my $sheet_name = 'combined_distance';
    my %option     = (
        chart_serial => 1,
        x_column     => 1,
        y_column     => 2,
        first_row    => 3,
        last_row     => 23,
        x_max_scale  => 20,
        x_title      => "Distance to indels (D1)",
        y_title      => "Nucleotide diversity",
        Height       => 200,
        Width        => 260,
        Top          => 14.25,
        Left         => 520,
    );

    #$excel_obj->draw_y( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- distance_non_slip
    #----------------------------#
    $sheet_name = 'distance_non_slip';

    $excel_obj->draw_y( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- distance_slip
    #----------------------------#
    $sheet_name = 'distance_slip';

    $excel_obj->draw_y( $sheet_name, \%option );

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

    #----------------------------#
    # worksheet -- align_te
    #----------------------------#
    my @te_levels = ( 1 .. 9 );

    foreach (@te_levels) {
        my $sheet_name = "align_te_$_";

        $excel_obj->draw_y( $sheet_name, \%option );
    }

    #----------------------------#
    # worksheet -- align_paralog
    #----------------------------#
    my @para_levels = ( 1 .. 9 );

    foreach (@para_levels) {
        my $sheet_name = "align_paralog_$_";

        $excel_obj->draw_y( $sheet_name, \%option );
    }

    #----------------------------#
    # worksheet -- distance_gc
    #----------------------------#
    $sheet_name = 'distance_gc';
    $option{y_title} = "GC proportion";

    #$excel_obj->draw_y( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- distance_gc_non_slip
    #----------------------------#
    $sheet_name = 'distance_gc_non_slip';

    $excel_obj->draw_y( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- distance_gc_slip
    #----------------------------#
    $sheet_name = 'distance_gc_slip';

    $excel_obj->draw_y( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- distance_dG
    #----------------------------#
    $sheet_name = 'distance_dG';
    $option{y_title} = "delta G";

    $excel_obj->draw_y( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- distance_dG_coding
    #----------------------------#
    $sheet_name = 'distance_dG_coding';

    $excel_obj->draw_y( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- distance_dG_non_coding
    #----------------------------#
    $sheet_name = 'distance_dG_non_coding';

    $excel_obj->draw_y( $sheet_name, \%option );

}

{

    #----------------------------#
    # worksheet -- combined_density
    #----------------------------#
    my $sheet_name = 'combined_density';
    my %option     = (
        chart_serial => 1,
        x_column     => 1,
        y_column     => 2,
        first_row    => 3,
        last_row     => 43,
        x_max_scale  => 20,
        x_title      => "Indel density (d2)",
        y_title      => "Nucleotide diversity",
        Height       => 200,
        Width        => 260,
        Top          => 14.25,
        Left         => 520,
    );

    #$excel_obj->draw_y( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- density_gc
    #----------------------------#
    $sheet_name = 'density_gc';
    $option{y_title} = "GC proportion";

    #$excel_obj->draw_y( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- density_dG
    #----------------------------#
    $sheet_name = 'density_dG';
    $option{y_title} = "delta G";

    $excel_obj->draw_y( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- density_dG_coding
    #----------------------------#
    $sheet_name = 'density_dG_coding';

    $excel_obj->draw_y( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- density_dG_non_coding
    #----------------------------#
    $sheet_name = 'density_dG_non_coding';

    $excel_obj->draw_y( $sheet_name, \%option );
}

{

    #----------------------------#
    # worksheet -- indel_size_group
    #----------------------------#
    # select worksheet by name
    my $sheet_name = 'indel_size_group';
    my @group_name = qw/1--5 6--10 11--50 51--300/;
    my %option     = (
        chart_serial   => 1,
        x_title        => "Distance to indels (d1)",
        y_title        => "Nucleotide diversity",
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

    #----------------------------#
    # worksheet -- indel_size_asymmetry
    #----------------------------#
    $sheet_name = 'indel_size_asymmetry';
    @group_name = (
        "left 1--10 & right 1--10",
        "left 1--10 & right 11--300",
        "left 11--300 & right 1--10",
        "left 11--300 & right 11--300",
    );
    $option{group_name} = \@group_name;

    $excel_obj->draw_dd( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- indel_extand_group
    #----------------------------#
    $sheet_name = 'indel_extand_group';
    @group_name = (
        "0--99",      "100--299", "300--499", "500--999",
        "1000--1999", "2000--99999",
    );
    $option{group_name} = \@group_name;

    $excel_obj->draw_dd( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- indel_extand_asymmetry
    #----------------------------#
    $sheet_name = 'indel_extand_asymmetry';
    @group_name = (
        "left 0--499 & right 0--499",
        "left 0--499 & right 500--99999",
        "left 500--99999 & right 0--499",
        "left 500--99999 & right 500--99999",
    );
    $option{group_name} = \@group_name;

    $excel_obj->draw_dd( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- indel_position_group
    #----------------------------#
    $sheet_name = 'indel_position_group';
    @group_name = (
        "coding & non-repeats",
        "coding & repeats",
        "non-coding & non-repeats",
        "non-coding & repeats",
    );
    $option{group_name} = \@group_name;

    $excel_obj->draw_dd( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- indel_coding_group
    #----------------------------#
    $sheet_name = 'indel_coding_group';
    @group_name = (
        "Deletion & non-coding",
        "Insertion & non-coding",
        "Deletion & coding",
        "Insertion & coding",
    );
    $option{group_name} = \@group_name;

    $excel_obj->draw_dd( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- indel_repeat_group
    #----------------------------#
    $sheet_name = 'indel_repeat_group';
    @group_name = (
        "Deletion & non-repeat",
        "Insertion & non-repeat",
        "Deletion & repeat",
        "Insertion & repeat",
    );
    $option{group_name} = \@group_name;

    $excel_obj->draw_dd( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- indel_slip_group
    #----------------------------#
    $sheet_name = 'indel_slip_group';
    @group_name = (
        "Deletion & non-slip",
        "Insertion & non-slip",
        "Deletion & slip",
        "Insertion & slip",
    );
    $option{group_name} = \@group_name;

    $excel_obj->draw_dd( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- indel_gc_group
    #----------------------------#
    $sheet_name = 'indel_gc_group';
    @group_name = (
        "Deletion & 0 <= gc < 0.3",
        "Insertion & 0 <= gc < 0.3",
        "Deletion & 0.3 <= gc < 0.5",
        "Insertion & 0.3 <= gc < 0.5",
        "Deletion & 0.5 <= gc <= 1",
        "Insertion & 0.5 <= gc <= 1",
    );
    $option{group_name} = \@group_name;

    $excel_obj->draw_dd( $sheet_name, \%option );
}

{

    #----------------------------#
    # worksheet -- snp_indel_ratio
    #----------------------------#
    my $sheet_name = 'snp_indel_ratio';
    my %option     = (
        chart_serial => 1,
        x_column     => 2,
        y_column     => 3,
        x_title      => "Nucleotide diversity",
        y_title      => "SNP/Indel ratio",
        Height       => 300,
        Width        => 390,
        Top          => 14.25,
        Left         => 620,
    );

    $excel_obj->draw_xy( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- snp_indel_ratio_2
    #----------------------------#
    $sheet_name = 'snp_indel_ratio_2';

    $excel_obj->draw_xy( $sheet_name, \%option );

}

{

    #----------------------------#
    # worksheet -- dd_group
    #----------------------------#
    # select worksheet by name
    my $sheet_name = 'dd_group';
    my @group_name
        = qw/200--399 400--799 800--1199 1200--1999 2000--2999 >=3000/;
    my %option = (
        chart_serial   => 1,
        x_title        => "Distance to indels (d1)",
        y_title        => "Nucleotide diversity",
        Height         => 300,
        Width          => 390,
        Top            => 14.25 * 38,
        Left           => 520,
        group_name     => \@group_name,
        section_top    => 2,
        section_end    => 36,
        section_length => 35,
        x_orientation  => 90,
    );

    $excel_obj->draw_dd( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- dd_group_gc
    #----------------------------#
    $sheet_name             = 'dd_group_gc';
    @group_name             = qw/200--599 600--1399 1400--2999 >=3000/;
    $option{group_name}     = \@group_name;
    $option{Top}            = 14.25 * 28;
    $option{section_top}    = 2;
    $option{section_end}    = 26;
    $option{section_length} = 25;
    $option{y_title}        = "GC proportion";

    $excel_obj->draw_dd( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- dd_group_dG
    #----------------------------#
    $sheet_name         = 'dd_group_dG';
    @group_name         = qw/200--399 400--1399 1400-2999 >=3000/;
    $option{group_name} = \@group_name;

    $excel_obj->draw_dd( $sheet_name, \%option );
}

{

    #----------------------------#
    # worksheet -- snp_base_change
    #----------------------------#
    my $sheet_name = 'snp_base_change';
    my @group_name = qw/T2Q Q2T Deletion Insertion/;
    my %option     = (
        chart_serial   => 1,
        y_scale_unit   => 2,
        x_title        => "Subsitution Type",
        y_title        => "Proportion of substitutions",
        Height         => 300,
        Width          => 390,
        Top            => 14.25 * 17,
        Left           => 360,
        group_name     => \@group_name,
        section_top    => 2,
        section_end    => 14,
        section_length => 13,
        x_orientation  => 90,
    );

    $excel_obj->draw_dd( $sheet_name, \%option );
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
        last_row      => 43,
        x_max_scale   => 40,
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
    $option{last_row}    = 23;
    $option{x_max_scale} = 20;
    $option{x_title}     = "Distance to indels (d1)";

    $excel_obj->draw_y( $sheet_name, \%option );
}

#----------------------------------------------------------#
# POST Processing
#----------------------------------------------------------#
# add time stamp to "summary" sheet
$excel_obj->time_stamp("basic") if $time_stamp;

# add an index sheet
$excel_obj->add_index_sheet if $add_index_sheet;

print "$outfile has been generated.\n";

$stopwatch->end_message;
exit;

__END__
    
=head1 NAME

    common_chart_factory.pl - Use Win32::OLE to automate Excel chart

=head1 SYNOPSIS

    common_chart_factory.pl [options]
     Options:
       --help           brief help message
       --man            full documentation
       --infile         input file name (full path)
       --outfile        output file name
       --jc             Jukes & Cantor correction
       --replace        replace text when charting
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
