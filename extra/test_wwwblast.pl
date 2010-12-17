#!/usr/bin/perl
use strict;
use warnings;

use WWW::Mechanize;
use Bio::SearchIO;
use YAML qw(Dump Load DumpFile LoadFile);

# Don't use http_proxy
our $saved_proxy = $ENV{'http_proxy'};
$ENV{'http_proxy'} = '';

$|++;

my $server         = 'http://202.119.43.222/blast/blast.html';
my $program        = 'blastn';
my $datalib        = 'Nipvs9311';
my $expact         = '0.0001';
my $alignment_view = '0';
my $filename       = "blast.out";
my $result_format  = 'blast';

my $sequence = q~
>test_seq
ATGCACTGTAGTAAAGGACTACTGCATCAGGCTGGGAGTCAAGATCTGCTTCGCTCGGTT
TCTCACCCTCAGAGCAATGGACAAGTCGAGAGGGCAAACGGCATAGTACTACAAGGAATC
AAGACCCGCGTCTACGACAGGCTCATGTCACATGACAAGAAGTGGGTCGAAGAACTTCCA
TCAGTACTATGGGCCGTGCGTACCACACCGACAACGTCTAACAAGGAGACACCCTTCTTC
CTCGTGTACGGCTCAGAGGCTATGCTCCCCACCAAGCTACGATACCAAAGTACACGAGCA
CATAAGTACTCCGATGAGAATCAGGAGGAGCAGCGAAATGACGATGTGAATCTACTCGAA
GAACATCGCGAACGAGTTGTCGTTCGAGCAGCCAGCTACCAACAGGCCCTTCACCGCTAC
TACGAGAAACGCATCCGAGCACGCACACTTTCGATCGGCGACTACGTCCTCCGACGGGTT
CAAAGTCAAGCAGGGCGAAACAAGCTCTCACCCAAATGGGAAGGACCGTACACGATCACA
CAAGTTCTGCAGCCAGGCGCATTCAAAATTGCAGACGGCGATGGTCGCGAGTTGGCAAAT
TCCTGGAACATTAATCAATTACGTAAATTCTATGTATAA
~;

my $seq_length = length $sequence;

# Fetch blast results section
{
    my $mech = WWW::Mechanize->new();

    print "Go to BLAST form page...\n";
    $mech->get($server);

    #print "Page content...\n";
    #print $mech->content();

    print "Submit query...\n";
    $mech->form_name('MainBlastForm');

    # Program:
    # <select name="PROGRAM">
    $mech->field( "PROGRAM", $program );

    # Expect:
    # <select name="DATALIB">
    $mech->field( "DATALIB", $datalib );

    # Database:
    # <select name="EXPECT">
    $mech->field( "EXPECT", $expact );

    # Filter:
    # <input name="FILTER">
    $mech->untick( "FILTER", "L" );

    # Alignment view:
    # <select name="ALIGNMENT_VIEW">
    $mech->field( "ALIGNMENT_VIEW", $alignment_view );

    # Sequence:
    # <textarea name="SEQUENCE">
    $mech->field( "SEQUENCE", $sequence );

    print "Waiting for server processing...\n";
    $mech->click();

    print "Save blast result to [$filename]...\n";
    my $blast_report = $mech->content( format => "text" );
    open FH, ">$filename";
    print FH $blast_report;
    close FH;

}

# Parse blast results section
{

    print "Reread in [$filename] file...\n";

    # Now parse blast file
    my $searchio = new Bio::SearchIO(
        -format => $result_format,
        -file   => $filename,
    );

    while ( my $result = $searchio->next_result() ) {
        while ( my $hit = $result->next_hit ) {

            # process the Bio::Search::Hit::HitI object
            my $hit_name     = $hit->name();
            my $desc         = $hit->description();
            my $length       = $hit->length;
            my $algorithm    = $hit->algorithm();
            my $score        = $hit->raw_score();
            my $significance = $hit->significance();
            my $rank = $hit->rank();    # the Nth hit for a specific query

            print Dump(
                {   hit_name     => $hit_name,
                    desc         => $desc,
                    length       => $length,
                    algorithm    => $algorithm,
                    score        => $score,
                    significance => $significance,
                    rank         => $rank,
                }
            );

            while ( my $hsp = $hit->next_hsp ) {

                # process the Bio::Search::HSP::HSPI object
                my $r_type = $hsp->algorithm;
                my $pvalue = $hsp->pvalue();
                my $evalue = $hsp->evalue();
                my $gaps   = $hsp->gaps( ['total'] );
                my $qseq   = $hsp->query_string;
                my $hseq   = $hsp->hit_string;
                my $len    = $hsp->length( ['total'] );
                my $rank   = $hsp->rank;
                my ( $hstart, $hend ) = $hsp->range( ['hit'] );

                print Dump(
                    {   r_type => $r_type,
                        pvalue => $pvalue,
                        evalue => $evalue,
                        gaps   => $gaps,
                        qseq   => $qseq,
                        hseq   => $hseq,
                        len    => $len,
                        rank   => $rank,
                        hstart => $hstart,
                        hend   => $hend,
                    }
                );
            }

        }
    }
}

print "Done.\n";

END {    # Restore $ENV{'http_proxy'}
    $ENV{'http_proxy'} = $saved_proxy;
}
