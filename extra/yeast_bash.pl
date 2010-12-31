#!/usr/bin/perl
use strict;
use warnings;

use Template;

{    # wgs
    my $tt   = Template->new;
    my @data = (
        { taxon => 226125, name => 'Spar',       coverage => '7x', },
        { taxon => 285006, name => 'RM11',       coverage => '10x', },
        { taxon => 307796, name => 'YJM789',     coverage => '10x', },
        { taxon => 574961, name => 'JAY291',     coverage => '162x', },
        { taxon => 538975, name => 'Sigma1278b', coverage => '45x', },
        { taxon => 545124, name => 'AWRI1631',   coverage => '7x', },
        { taxon => 643680, name => 'EC1118',     coverage => 'unknown', },
        { taxon => 538976, name => 'YPS163',     coverage => '2.8x', },
        { taxon => 538975, name => 'M22',        coverage => '2.6x', },
    );

    my $text = <<'EOF';
#!/bin/bash
# cd ~/data/alignment/yeast

[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
gunzip -c [% item.name %]/*.gz > [% item.name %]/[% item.name %].fasta
perl -p -i -e '/>/ and s/\>gi\|(\d+).*/\>gi_$1/' [% item.name %]/[% item.name %].fasta
RepeatMasker [% item.name %]/*.fasta -species fungi -xsmall -s --parallel 8
mv [% item.name %]/[% item.name %].fasta.masked [% item.name %]/[% item.name %].fa
perl /home/wangq/Scripts/blastz/bz.pl -dt /home/wangq/data/alignment/yeast/S288C_58 -dq /home/wangq/data/alignment/yeast/[% item.name %] -dl /home/wangq/data/alignment/yeast/S288Cvs[% item.name %]_58 -s set02 -p 6 
perl /home/wangq/Scripts/alignDB/extra/two_way_batch.pl -d S288Cvs[% item.name %] -e yeast_58 -t="4932,S288C" -q "[% item.taxon %],[% item.name %]" -a /home/wangq/data/alignment/yeast/S288Cvs[% item.name %]_58 -at 10000 -st 1000000 --parallel 4 --run all

[% END -%]
EOF

    $tt->process( \$text, { data => \@data }, "auto_wgs.sh" )
        or die Template->error;
}

{    # sgrp
    my $tt   = Template->new;
    my @data = (
        { taxon => 900001, name => 'Y55',       coverage => '4.1x', },
        { taxon => 580239, name => 'SK1',       coverage => '3.8x', },
        { taxon => 580240, name => 'W303',      coverage => '3.7x', },
        { taxon => 900003, name => 'DBVPG6765', coverage => '3.5x', },
    );

    my $text = <<'EOF';
#!/bin/bash
# cd ~/data/alignment/yeast

[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
RepeatMasker [% item.name %]/*.fasta -species fungi -xsmall -s --parallel 8
mv [% item.name %]/[% item.name %].fasta.masked [% item.name %]/[% item.name %].fa
perl /home/wangq/Scripts/blastz/bz.pl -dt /home/wangq/data/alignment/yeast/S288C_58 -dq /home/wangq/data/alignment/yeast/[% item.name %] -dl /home/wangq/data/alignment/yeast/S288Cvs[% item.name %]_58 -s set02 -p 6 
perl /home/wangq/Scripts/alignDB/extra/two_way_batch.pl -d S288Cvs[% item.name %] -e yeast_58 -t="4932,S288C" -q "[% item.taxon %],[% item.name %]" -a /home/wangq/data/alignment/yeast/S288Cvs[% item.name %]_58 -at 10000 -st 1000000 --parallel 4 --run all

[% END -%]
EOF

    $tt->process( \$text, { data => \@data }, "auto_sgrp.sh" )
        or die Template->error;
}

{    # join dbs
    my $tt   = Template->new;
    my @data = (
        {   goal_db  => 'S288CvsThree',
            dbs      => 'S288CvsSpar,S288CvsRM11,S288CvsYJM789',
            outgroup => '0query',
            target   => '0target',
            queries  => '1query,2query',
        },
        {   goal_db => 'S288CvsSix',
            dbs     => 'S288CvsSpar,S288CvsRM11,S288CvsYJM789'
                . ',S288CvsJAY291,S288CvsSigma1278b,S288CvsAWRI1631',
            outgroup => '0query',
            target   => '0target',
            queries  => '1query,2query,3query,4query,5query',
        },
        {   goal_db => 'S288CvsFourteen',
            dbs     => 'S288CvsSpar,S288CvsRM11,S288CvsYJM789'
                . ',S288CvsJAY291,S288CvsSigma1278b,S288CvsAWRI1631'
                . ',S288CvsEC1118,S288CvsYPS163,S288CvsM22'
                . ',S288CvsY55,S288CvsSK1,S288CvsW303,S288CvsDBVPG6765',
            outgroup => '0query',
            target   => '0target',
            queries  => '1query,2query,3query,4query,5query'
                . ',6query,7query,8query,9query,10query,11query,12query',
        },
    );

    my $text = <<'EOF';
#!/bin/bash
# cd ~/data/alignment/yeast

[% FOREACH item IN data -%]
# [% item.goal_db %]
perl /home/wangq/Scripts/alignDB/extra/join_dbs.pl --dbs [% item.dbs %] \
    --goal_db [% item.goal_db _ "_10k" %] --outgroup [% item.outgroup %] --target [% item.target %] \
    --queries [% item.queries %] \
    --no_insert=1 --trimmed_fasta=1 --length 10000

perl /home/wangq/Scripts/alignDB/extra/join_dbs.pl --dbs [% item.dbs %] \
    --goal_db [% item.goal_db _ "_5k" %] --outgroup [% item.outgroup %] --target [% item.target %] \
    --queries [% item.queries %] \
    --no_insert=1 --trimmed_fasta=1 --length 5000

[% END -%]
EOF

    $tt->process( \$text, { data => \@data }, "auto_joins.sh" )
        or die Template->error;
}
