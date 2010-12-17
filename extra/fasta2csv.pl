#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use List::Util qw(first max maxstr min minstr reduce shuffle sum);
use List::MoreUtils qw(any all);

#----------------------------------------------------------#
# GetOpt Section
#----------------------------------------------------------#
my $infile        = "../116_112.fasta";
my $wrap          = 60;
my $spacing       = 0;
my $outfile       = "";
my $has_reference = 1;
my $filter        = "ingroup";            # outgroup or ingroup

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'      => \$help,
    'man'         => \$man,
    'infile=s'    => \$infile,
    'wrap=i'      => \$wrap,
    'spacing=i'   => \$spacing,
    'outfile=s'   => \$outfile,
    'reference=s' => \$has_reference,
    'filter=s'    => \$filter,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

unless ($outfile) {
    $outfile = $infile;
    if ( $outfile =~ /\./ ) {
        $outfile =~ s/\.\w*?$/\.csv/;
    }
    else {
        $outfile .= ".csv";
    }
}

#----------------------------------------------------------#
# START, read in seqs
#----------------------------------------------------------#
# Read out file
open my $fasta_fh, '<', $infile
    or die "Cannot open Pearson file $infile";
my @file_contents = <$fasta_fh>;
close $fasta_fh;

my @seq_names;
my %seqs;
foreach my $current_line (@file_contents) {
    if ( $current_line =~ /^\>([\w-])+/ ) {
        $current_line =~ s/\>//;
        chomp $current_line;
        push @seq_names, $current_line;
        $seqs{$current_line} = '';
    }
    elsif ( $current_line =~ /^[\w-]+/ ) {
        $current_line =~ s/[^\w-]//g;
        chomp $current_line;
        my $current_seq_name = $seq_names[-1];
        $seqs{$current_seq_name} .= $current_line;
    }
    else { # Blank line, do nothing
    }
}

# seq names
my $ref_name;
my @total_names;
if ($has_reference) {
    @total_names = @seq_names;
    $ref_name    = shift @seq_names;
}
else {
    @total_names = @seq_names;
}

# seq number
my $seq_number   = scalar @seq_names;
my $total_number = scalar @total_names;

# seq length
my @total_lengths;
foreach (@total_names) {
    push @total_lengths, length $seqs{$_};
}
unless ( all { $_ == $total_lengths[0] } @total_lengths ) {
    die "Seq length error.\n";
}
else {
    print "Length: $total_lengths[0] bp\n";
}

#----------------------------------------------------------#
# Collect SNP sites (indel and Ns are treated as SNPs)
#----------------------------------------------------------#
my $result_matrix = [ [] ];    # csv matrix
my $section       = 0;
my $snp_count     = 0;

# add seq names to lines
foreach ( 0 .. $total_number ) {
    push @{ $result_matrix->[$_] }, ( " ", @total_names )[$_];
}
for ( my $i = 1; $i <= $total_lengths[0]; $i++ ) {
    my @cur_column;
    my @cur_in_column;
    foreach (@total_names) {
        my $base = substr( $seqs{$_}, $i - 1, 1 );
        push @cur_column, $base;
    }

    # if ingroup, reference is not count as SNP sites
    if ( $filter eq "ingroup" ) {
        @cur_in_column = @cur_column;
        shift @cur_in_column;
    }
    else {
        @cur_in_column = @cur_column;
    }

    if ( minstr(@cur_in_column) ne maxstr(@cur_in_column) ) {
        unshift @cur_column, $i;
        $snp_count++;
        foreach ( 0 .. $total_number ) {
            push @{   $result_matrix->[ $_
                    + $section * ( $total_number + 1 + $spacing ) ] },
                $cur_column[$_];
        }

        # add seq names to wrapped lines
        if ( ( $snp_count % $wrap ) == 0 ) {
            $section++;
            print "Section: $section\n";
            foreach ( 0 .. $total_number ) {
                push @{   $result_matrix->[ $_
                        + $section * ( $total_number + 1 + $spacing ) ] },
                    ( " ", @total_names )[$_];
            }
        }
    }
}

#----------------------------------------------------------#
# Output CSV
#----------------------------------------------------------#
open my $csv_fh, '>', $outfile
    or die "Cannot write CSV file $outfile";
foreach (@$result_matrix) {
    if ( ref($_) eq 'ARRAY' ) {
        print {$csv_fh} join( ",", (@$_) ), "\n";
    }
    else {
        print {$csv_fh} "\n";
    }
}
close $csv_fh;

__END__

=head1 NAME

    fasta2csv.pl - Extract polymorphic sites from fasta alignment file

=head1 SYNOPSIS

    fasta2csv.pl [options]
     Options:
       --help            brief help message
       --man             full documentation
       --infile          input file name (full path)
       --wrap            wrap number 
       --spacing         wrapped line spacing
       --outfile         output file name
       --reference       has reference?
       --filter          ingroup or outgroup sites
       

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
