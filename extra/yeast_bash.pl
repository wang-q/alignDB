#!/usr/bin/perl
use strict;
use warnings;

use Template;

my $tt = Template->new;

# wgs
my @strains_wgs = (
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

my $tt_wgs = <<'EOF';
#!/bin/bash
# cd ~/data/alignment/yeast

[% FOREACH strain IN strains -%]
# [% strain.name %] [% strain.coverage %]
gunzip -c [% strain.name %]/*.gz > [% strain.name %]/[% strain.name %].fasta
perl -p -i -e '/>/ and s/\>gi\|(\d+).*/\>gi_$1/' [% strain.name %]/[% strain.name %].fasta
RepeatMasker [% strain.name %]/*.fasta -species fungi -xsmall -s --parallel 8
mv [% strain.name %]/[% strain.name %].fasta.masked [% strain.name %]/[% strain.name %].fa
perl /home/wangq/Scripts/blastz/bz.pl -dt /home/wangq/data/alignment/yeast/S288C_58 -dq /home/wangq/data/alignment/yeast/[% strain.name %] -dl /home/wangq/data/alignment/yeast/S288Cvs[% strain.name %]_58 -s set02 -p 6 
perl /home/wangq/Scripts/alignDB/extra/two_way_batch.pl -d S288Cvs[% strain.name %] -e yeast_58 -t="4932,S288C" -q "[% strain.taxon %],[% strain.name %]" -a /home/wangq/data/alignment/yeast/S288Cvs[% strain.name %]_58 -at 10000 -st 1000000 --parallel 4 --run all

[% END -%]
EOF

$tt->process( \$tt_wgs, { strains => \@strains_wgs }, "auto_wgs.sh" )
    or die Template->error;

# sgrp
my @strains_sgrp = (
    { taxon => 900001, name => 'Y55',       coverage => '4.1x', },
    { taxon => 580239, name => 'SK1',       coverage => '3.8x', },
    { taxon => 580240, name => 'W303',      coverage => '3.7x', },
    { taxon => 900003, name => 'DBVPG6765', coverage => '3.5x', },
);

my $tt_sgrp = <<'EOF';
#!/bin/bash
# cd ~/data/alignment/yeast

[% FOREACH strain IN strains -%]
# [% strain.name %] [% strain.coverage %]
RepeatMasker [% strain.name %]/*.fasta -species fungi -xsmall -s --parallel 8
mv [% strain.name %]/[% strain.name %].fasta.masked [% strain.name %]/[% strain.name %].fa
perl /home/wangq/Scripts/blastz/bz.pl -dt /home/wangq/data/alignment/yeast/S288C_58 -dq /home/wangq/data/alignment/yeast/[% strain.name %] -dl /home/wangq/data/alignment/yeast/S288Cvs[% strain.name %]_58 -s set02 -p 6 
perl /home/wangq/Scripts/alignDB/extra/two_way_batch.pl -d S288Cvs[% strain.name %] -e yeast_58 -t="4932,S288C" -q "[% strain.taxon %],[% strain.name %]" -a /home/wangq/data/alignment/yeast/S288Cvs[% strain.name %]_58 -at 10000 -st 1000000 --parallel 4 --run all

[% END -%]
EOF

$tt->process( \$tt_wgs, { strains => \@strains_sgrp }, "auto_sgrp.sh" )
    or die Template->error;
