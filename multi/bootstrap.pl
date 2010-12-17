#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use Spreadsheet::ParseExcel;

use FindBin;
use lib "$FindBin::Bin/../lib";
use AlignDB::Util qw(:all);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

my $in_file;
my $rep        = 10000;             # times of bootstrap simulations
my $freq       = 1;                 # indel frequence, 1/3 or 2/3
my $sheet_name = "one_third_all";

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'  => \$help,
    'man'     => \$man,
    'in=s'    => \$in_file,
    'rep=i'   => \$rep,
    'freq=s'  => \$freq,
    'sheet=s' => \$sheet_name,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Start
#----------------------------------------------------------#
my $dispatch_func = {
    1 => \&mutation_rate_1,
    2 => \&mutation_rate_2,
    all => \&mutation_rate_all,
};
my $dispatch_range = {
    1 => \&range_1,
    2 => \&range_2,
    all => \&range_all,
};

my $excel = Spreadsheet::ParseExcel::Workbook->Parse($in_file);

foreach my $sheet ( @{ $excel->{Worksheet} } ) {
    next unless $sheet->{Name} eq $sheet_name;

    print "Worksheet $sheet->{Name}\n";

    my @data;
    foreach my $row ( $sheet->{MinRow} + 1 .. $sheet->{MaxRow} ) {
        my %info = (
            id  => $sheet->{Cells}[$row][0]{Val},
            S_i => $sheet->{Cells}[$row][5]{Val},
            S_n => $sheet->{Cells}[$row][6]{Val},
        );
        push @data, \%info;
    }

    my $ori_Ni_Nn  = calc_ratio( \@data );
    my $ori_uhet_u = $dispatch_func->{$freq}->($ori_Ni_Nn);

    print "Original Ni/Nn is ",  $ori_Ni_Nn,  "\n";
    print "Original uhet/u is ", $ori_uhet_u, "\n";

    my ( @Ni_Nn, @uhet_u );
    for ( 1 .. $rep ) {
        my $bs_data = sampling_with_replacement( \@data );

        my $bs_Ni_Nn = calc_ratio($bs_data);
        redo unless $dispatch_range->{$freq}->($bs_Ni_Nn);

        my $bs_uhet_u = $dispatch_func->{$freq}->($bs_Ni_Nn);

        #push @Ni_Nn,  $bs_Ni_Nn;
        push @uhet_u, $bs_uhet_u;
    }

    #print "Stats of bootstrap Ni/Nn is\n",  Dump stat_result( \@Ni_Nn );
    print "Stats of bootstrap uhet/u is\n", Dump stat_result( \@uhet_u );

    print "\n" x 2;
}

exit;

sub calc_ratio {
    my $data = shift;

    my ( $sum_S_i, $sum_S_n );
    for (@$data) {
        $sum_S_i += $_->{S_i};
        $sum_S_n += $_->{S_n};
    }

    my $ratio = $sum_S_i / $sum_S_n;

    return $ratio;
}

sub mutation_rate_1 {
    my $ratio = shift;
    return ( 17 - 21 * $ratio ) / ( 3 * $ratio - 7 );
}

sub mutation_rate_2 {
    my $ratio = shift;
    return ( 36 - 44 * $ratio ) / ( 12 * $ratio - 20 );
}

sub mutation_rate_all {
    my $ratio = shift;
    return ( 35 - 43 * $ratio ) / ( 9 * $ratio - 17 );
}

sub range_1 {
    my $ratio = shift;
    if ( $ratio >= 17/21 and $ratio < 7 / 3 ) {
        return 1;
    }
    else {
        return 0;
    }
}

sub range_2 {
    my $ratio = shift;
    if ( $ratio >= 9/11 and $ratio < 5 / 3 ) {
        return 1;
    }
    else {
        return 0;
    }
}

sub range_all {
    my $ratio = shift;
    if ( $ratio >= 35/43 and $ratio < 17 / 9 ) {
        return 1;
    }
    else {
        return 0;
    }
}

__END__

=head1 NAME

    bootstrap.pl - bootstrap the indel-induced mutation rate increase

=head1 SYNOPSIS

    apply_sql.pl [options]
     Options:
       --help            brief help message
       --man             full documentation
       --file            csv file name
       --rep             times of bootstrap simulations

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

