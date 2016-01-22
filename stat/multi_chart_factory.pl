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

if (! -e $infile) {
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
my @sheet_names = @{ $excel_obj->sheet_names };
{

    #----------------------------#
    # worksheet -- indel_type_gc_10
    #----------------------------#
    my $sheet_name = 'indel_type_gc_10';
    my @group_name = qw{Insertion Deletion};
    my %option     = (
        chart_serial   => 1,
        x_title        => "Indel length",
        y_title        => "GC proportion",
        Height         => 300,
        Width          => 390,
        Top            => 14.25 * 17,
        Left           => 360,
        group_name     => \@group_name,
        section_top    => 2,
        section_end    => 12,
        section_length => 11,
        x_orientation  => 0,
    );

    $excel_obj->draw_dd( $sheet_name, \%option );

    #----------------------------#
    # worksheet -- indel_type_gc_100
    #----------------------------#
    $sheet_name = 'indel_type_gc_100';

    $excel_obj->draw_dd( $sheet_name, \%option );
}

#----------------------------------------------------------#
# POST Processing
#----------------------------------------------------------#

print "$outfile has been generated.\n";

$stopwatch->end_message;
exit;

__END__
    
=head1 NAME

    multi_chart_factory.pl - Use Win32::OLE to automate Excel chart

=head1 SYNOPSIS

    multi_chart_factory.pl [options]
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
