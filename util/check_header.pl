#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Getopt::Long;
use Pod::Usage;
use YAML qw(Dump Load DumpFile LoadFile);

use File::Spec;
use File::Basename;
use IO::Zlib;

use AlignDB::Stopwatch;
use AlignDB::Util qw(:all);

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#
my $in_file;    # be checked file, normal multi fasta or blocked fasta
my $genome;     # multi fasta containing the genome
my $name;       # which species to be checked (should be same as $genome)
                # omit this will check all sequences

my $samtools = 'samtools';    # samtools exec file
my $detail;                   # log error sequences
my $gzip;                     # open .gz

my $man  = 0;
my $help = 0;

GetOptions(
    'help|?'       => \$help,
    'man'          => \$man,
    'i|in_file=s'  => \$in_file,
    'g|genome=s'   => \$genome,
    'n|name=s'     => \$name,
    's|samtools=s' => \$samtools,
    'detail'       => \$detail,
    'gzip'         => \$gzip,
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

#----------------------------------------------------------#
# Search for all files
#----------------------------------------------------------#
my $stopwatch = AlignDB::Stopwatch->new;

my $in_fh;
if ( !$gzip ) {
    open $in_fh, '<', $in_file;
}
else {
    $in_fh = IO::Zlib->new( $in_file, "rb" );
}

my $log_fh;
if ($detail) {
    open $log_fh, '>', basename($in_file) . '.log.txt';
}

{
    my $header;
    my $content = '';
    my $count   = 0;
    while ( my $line = <$in_fh> ) {
        chomp $line;

        if ( $line =~ /^\>[\w:-]+/ ) {

            # the first sequence is ready
            if ( defined $header ) {
                $count = check_seq( $header, $content, $count );
            }

            # prepare to accept next sequence
            $line =~ s/^\>//;
            $header = $line;

            # clean previous sequence
            $content = '';
        }
        elsif ( $line =~ /^[\w-]+/ ) {
            $line =~ s/[^\w]//g;    # including '-'s
            $line = uc $line;
            $content .= $line;
        }
        else {                      # Blank line, do nothing
        }
    }
    
    # for last sequece
    check_seq( $header, $content, $count );
}

if ( !$gzip ) {
    close $in_fh;
}
else {
    $in_fh->close;
}

if ($detail) {
    close $log_fh;
}

$stopwatch->block_message( "All files have been processed.", "duration" );

exit;

sub get_seq_faidx {
    my $samtools = shift;
    my $genome   = shift;
    my $location = shift;

    my $cmd = sprintf "%s faidx %s %s", $samtools, $genome, $location;
    open my $fh_pipe, '-|', $cmd;

    my $seq;
    while ( my $line = <$fh_pipe> ) {
        chomp $line;
        if ( $line =~ /^[\w-]+/ ) {
            $seq .= $line;
        }
    }
    close($fh_pipe);

    return $seq;
}

sub check_seq {
    my $header = shift;
    my $seq    = shift;
    my $count  = shift;

    $count++;
    my $info = decode_header($header);
    if ( $name and $name ne $info->{name} ) {
        next;
    }

    if ( $info->{chr_strand} eq '-' ) {
        $seq = revcom($seq);
    }

    my $location;
    if ( $info->{chr_end} ) {
        $location = sprintf "%s:%s-%s", $info->{chr_name},
            $info->{chr_start}, $info->{chr_end};
    }
    else {
        $location = sprintf "%s:%s", $info->{chr_name}, $info->{chr_start};
    }
    my $seq_genome = get_seq_faidx( $samtools, $genome, $location );

    if ( $seq eq $seq_genome ) {
        printf "OK\t%s\t%s\n", $count, $header;
    }
    else {
        printf "FAILED\t%s\t%s\n", $count, $header;
        if ($detail) {
            print {$log_fh} ">$header\n";
            print {$log_fh} "$seq\n";
            print {$log_fh} ">$location\n";
            print {$log_fh} "$seq_genome\n";
            print {$log_fh} "\n";
        }
    }

    return $count;
}

__END__

=head1 NAME

    check_header.pl - check genome location in (blocked) fasta headers

=head1 SYNOPSIS
    perl check_header.pl --in I.net.axt.fas -g S288c.fasta --detail

    check_header.pl [options]
      Options:
        --help              brief help message
        --man               full documentation
        --in                axt file's location
        --genome            one multi fasta file contains genome sequences
                            extract sub-sequence by samtools faidx
        --detail            write a fasta file report error sequences
        --gzip              input file is gzipped

=cut

find ~/data/alignment/yeast_genome/S288c/ -name "*.fasta" \
    | sort \
    | xargs cat > ~/data/alignment/yeast_genome/S288c.fasta

perl check_header.pl --in ~/data/alignment/self/yeast_new/S288cvsselfalign_fasta/I.net.axt.gz.fas \
    -g ~/data/alignment/yeast_genome/S288c.fasta \
    --detail
