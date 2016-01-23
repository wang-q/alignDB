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
my $Config = Config::Tiny->new();
$Config = Config::Tiny->read("$FindBin::Bin/../alignDB.ini");

# stat parameter
my $jc_correction   = $Config->{stat}->{jc_correction};
my $time_stamp      = $Config->{stat}->{time_stamp};
my $add_index_sheet = $Config->{stat}->{add_index_sheet};

# in three-lineage test, we don't want jc_correction
$jc_correction = 0;

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
my $stopwatch = AlignDB::Stopwatch->new();
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
    # worksheet -- distance
    #----------------------------#
    my $sheet_name = 'distance';
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

    #----------------------------#
    # worksheet -- distance_window0eq0
    #----------------------------#
    $sheet_name            = 'distance_window0eq0';
    $option{chart_serial}  = 1;
    $option{last_row}      = 7;
    $option{x_max_scale}   = 5;
    $option{y_column}      = 3;
    $option{y_last_column} = 4;
    $option{y_title}       = "Nucleotide diversity";
    $option{y_scale_unit}  = 0.001;
    $option{Top}           = 12.75;

    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 2
    $option{chart_serial}++;
    $option{y_column}      = 7;
    $option{y_last_column} = 7;
    $option{y_title}       = "Di/Dn";
    $option{y_scale_unit}  = 0.1;
    $option{Top} += $option{Height} + 12.75;

    $excel_obj->draw_y( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- distance_indel_bottom90
    #----------------------------#
    $sheet_name            = 'distance_indel_bottom90';
    $option{chart_serial}  = 1;
    $option{last_row}      = 7;
    $option{x_max_scale}   = 5;
    $option{y_column}      = 3;
    $option{y_last_column} = 4;
    $option{y_title}       = "Nucleotide diversity";
    $option{y_scale_unit}  = 0.001;
    $option{Top}           = 12.75;

    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 2
    $option{chart_serial}++;
    $option{y_column}      = 7;
    $option{y_last_column} = 7;
    $option{y_title}       = "Di/Dn";
    $option{y_scale_unit}  = 0.1;
    $option{Top} += $option{Height} + 12.75;

    $excel_obj->draw_y( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- combined_distance
    #----------------------------#
    $sheet_name            = 'combined_distance';
    $option{chart_serial}  = 1;
    $option{last_row}      = 32;
    $option{x_max_scale}   = 30;
    $option{y_column}      = 3;
    $option{y_last_column} = 4;
    $option{y_title}       = "Nucleotide diversity";
    $option{y_scale_unit}  = 0.001;
    $option{Top}           = 12.75;

    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 2
    $option{chart_serial}++;
    $option{y_column}      = 7;
    $option{y_last_column} = 7;
    $option{y_title}       = "Di/Dn";
    $option{y_scale_unit}  = 0.1;
    $option{Top} += $option{Height} + 12.75;

    $excel_obj->draw_y( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- distance_coding
    #----------------------------#
    $sheet_name            = 'distance_coding';
    $option{chart_serial}  = 1;
    $option{last_row}      = 32;
    $option{x_max_scale}   = 30;
    $option{y_column}      = 3;
    $option{y_last_column} = 4;
    $option{y_title}       = "Nucleotide diversity";
    $option{y_scale_unit}  = 0.001;
    $option{Top}           = 12.75;

    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 2
    $option{chart_serial}++;
    $option{y_column}      = 7;
    $option{y_last_column} = 7;
    $option{y_title}       = "Di/Dn";
    $option{y_scale_unit}  = 0.1;
    $option{Top} += $option{Height} + 12.75;

    $excel_obj->draw_y( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- distance_non_coding
    #----------------------------#
    $sheet_name            = 'distance_non_coding';
    $option{chart_serial}  = 1;
    $option{last_row}      = 32;
    $option{x_max_scale}   = 30;
    $option{y_column}      = 3;
    $option{y_last_column} = 4;
    $option{y_title}       = "Nucleotide diversity";
    $option{y_scale_unit}  = 0.001;
    $option{Top}           = 12.75;

    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 2
    $option{chart_serial}++;
    $option{y_column}      = 7;
    $option{y_last_column} = 7;
    $option{y_title}       = "Di/Dn";
    $option{y_scale_unit}  = 0.1;
    $option{Top} += $option{Height} + 12.75;

    $excel_obj->draw_y( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- snp_indel
    #----------------------------#
    $sheet_name            = 'snp_indel';
    $option{chart_serial}  = 1;
    $option{last_row}      = 12;
    $option{x_max_scale}   = 10;
    $option{y_column}      = 2;
    $option{y_last_column} = 2;
    $option{y_title}       = "Indel";
    $option{x_title}       = "Distance to snp";
    $option{y_scale_unit}  = 0.01;
    $option{Top}           = 12.75;

    $excel_obj->draw_y( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- ds_dns
    #----------------------------#
    $sheet_name            = 'ds_dns';
    $option{chart_serial}  = 1;
    $option{last_row}      = 31;
    $option{x_max_scale}   = 30;
    $option{chart_serial}  = 1;
    $option{y_column}      = 2;
    $option{y_last_column} = 3;
    $option{x_title}       = "Distance to interval center";
    $option{y_title}       = "Nucleotide diversity";
    $option{Top}           = 12.75;

    $excel_obj->draw_y( $sheet_name, \%option );

    # chart 2
    $option{chart_serial}++;
    $option{y_column}      = 5;
    $option{y_last_column} = 5;
    $option{y_title}       = "Ds/Dns";
    $option{y_scale_unit}  = 0.1;
    $option{Top} += $option{Height} + 12.75;

    $excel_obj->draw_y( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- distance_cpg
    #----------------------------#
    $sheet_name            = 'distance_cpg';
    $option{chart_serial}  = 1;
    $option{first_row}     = 3;
    $option{last_row}      = 33;
    $option{x_max_scale}   = 30;
    $option{chart_serial}  = 1;
    $option{y_column}      = 2;
    $option{y_last_column} = 2;
    $option{x_title}       = "Distance to indel";
    $option{y_title}       = "CpG/100bp";
    $option{Top}           = 12.75;
    $option{Left}          = 220;

    $excel_obj->draw_y( $sheet_name, \%option );

}

{

    #----------------------------#
    # worksheet -- Dxr
    #----------------------------#
    my @sheets = qw{
        distance_dir_dnr
        distance_dtr_dqr
        target_dtr_dqr
        query_dtr_dqr
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
    # worksheet -- snp_distance
    #----------------------------#
    foreach (
        qw/all -1 0-5 6-10 11-20 21-30 31-40 41-50 51-999/,
        qw/21-999 31-999 41-999/,
        qw/0-0 1-1 2-2 3-3 4-4 5-5 6-6 7-7 8-8 9-9 10-10/
        )
    {
        my $sheet_name = 'snp_distance_' . $_;
        my %option     = (
            chart_serial  => 1,
            x_column      => 1,
            y_column      => 3,
            y_last_column => 4,
            first_row     => 2,
            last_row      => 12,
            x_max_scale   => 10,
            y_scale_unit  => 0.001,
            x_title       => "Distance to snp",
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
        $option{y_title}       = "Ds/Dns";
        $option{y_scale_unit}  = 0.1;
        $option{Top} += $option{Height} + 12.75;
        $excel_obj->draw_y( $sheet_name, \%option );
    }

}

{

    #----------------------------#
    # worksheet -- snp_LR_distance
    #----------------------------#
    foreach (
        qw/R-0-0 R-1-1 R-2-2 R-3-3 R-4-4 R-5-5/,
        qw/L-0-0 L-1-1 L-2-2 L-3-3 L-4-4 L-5-5/,
        )
    {
        my $sheet_name = 'snp_distance_' . $_;
        my %option     = (
            chart_serial  => 1,
            x_column      => 1,
            y_column      => 2,
            y_last_column => 2,
            first_row     => 1,
            last_row      => 25,
            y_scale_unit  => 0.001,
            x_title       => "Distance to snp",
            y_title       => "Nucleotide diversity",
            Height        => 283.7,
            Width         => 453.9,
            Top           => 12.75,
            Left          => 440,
        );
        $excel_obj->draw_LineMarkers( $sheet_name, \%option );

        ## chart 2
        #$option{chart_serial}++;
        #$option{y_column}      = 7;
        #$option{y_last_column} = 7;
        #$option{y_title}       = "Ds/Dns";
        #$option{y_scale_unit}  = 0.1;
        #$option{Top} += $option{Height} + 12.75;
        #$excel_obj->draw_y($sheet_name, \%option);
    }

}

{

    #----------------------------#
    # worksheet -- indel_distance_
    #----------------------------#
    foreach (qw/non-slip slip insertion deletion target query/) {
        my $sheet_name = 'indel_distance_' . $_;
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

{

    #----------------------------#
    # worksheet -- indel_slip_di_group
    #----------------------------#
    # select worksheet by name
    my $sheet_name = 'indel_slip_di_group';
    my @group_name = qw/0--0 1--1/;
    my %option     = (
        chart_serial   => 1,
        y_scale_unit   => 0.001,
        x_title        => "Distance to indels (d1)",
        y_title        => "Nucleotide diversity",
        Height         => 283.7,
        Width          => 453.9,
        Top            => 12.75 * 17,
        Left           => 360,
        group_name     => \@group_name,
        section_top    => 2,
        section_end    => 15,
        section_length => 14,
        x_orientation  => 0,
    );

    $excel_obj->draw_dd( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- indel_slip_dn_group
    #----------------------------#
    # select worksheet by name
    $sheet_name = 'indel_slip_dn_group';
    @group_name = qw/0--0 1--1/;

    $excel_obj->draw_dd( $sheet_name, \%option );
}

{

    #----------------------------#
    # worksheet -- snp_base_change
    #----------------------------#
    my $sheet_name = 'snp_base_change';
    my @group_name = qw/Target Query/;
    my %option     = (
        chart_serial   => 1,
        y_scale_unit   => 2,
        x_title        => "Subsitution Type",
        y_title        => "Proportion of substitutions",
        Height         => 283.7,
        Width          => 453.9,
        Top            => 12.75 * 17,
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
    # worksheet -- distance_snp
    #----------------------------#
    my $sheet_name = 'distance_snp';
    my %option     = (
        chart_serial  => 1,
        x_column      => 1,
        y_column      => 2,
        y_last_column => 13,
        first_row     => 1,
        last_row      => 13,
        x_min_scale   => -1,
        x_max_scale   => 10,
        y_scale_unit  => 0.05,
        x_title       => "Distance to indels (d1)",
        y_title       => "Proportion of substitutions",
        x_ori         => 0,
        Height        => 283.7,
        Width         => 453.9,
        Top           => 12.75 * 2,
        Left          => 200,
    );

    $excel_obj->draw_LineMarkers( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- distance_snp_non_cpg
    #----------------------------#
    $sheet_name = 'distance_snp_non_cpg';

    $excel_obj->draw_LineMarkers( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- distance_tri_trv
    #----------------------------#
    $sheet_name = 'distance_tri_trv';
    $option{y_last_column} = 3;

    $excel_obj->draw_LineMarkers( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- distance_tri_trv_non_cpg
    #----------------------------#
    $sheet_name = 'distance_tri_trv_non_cpg';
    $excel_obj->draw_LineMarkers( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- distance_tri_trv_coding
    #----------------------------#
    $sheet_name = 'distance_tri_trv_coding';

    $excel_obj->draw_LineMarkers( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- distance_tri_trv_non_coding
    #----------------------------#
    $sheet_name = 'distance_tri_trv_non_coding';

    $excel_obj->draw_LineMarkers( $sheet_name, \%option );
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

three_chart_factory.pl - Use Win32::OLE to automate Excel chart

=head1 SYNOPSIS

    three_chart_factory.pl [options]
     Options:
       --help            brief help message
       --man             full documentation
       --infile          input file name (full path)
       --outfile         output file name

=cut
