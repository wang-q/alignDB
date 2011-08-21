#!/usr/bin/perl
use strict;
use warnings;

use Template;

my $data_dir = shift @ARGV;
$data_dir ||= "/home/wangq/data/alignment/yeast";

my $pl_dir = shift @ARGV;
$pl_dir ||= "/home/wangq/Scripts";

{    # wgs
    my $tt   = Template->new;
    my @data = (
        { taxon => 226125, name => 'Spar',       coverage => '7x', },
        { taxon => 285006, name => 'RM11',       coverage => '10x', },
        { taxon => 307796, name => 'YJM789',     coverage => '10x', },
        { taxon => 574961, name => 'JAY291',     coverage => '162x', },
        { taxon => 538975, name => 'Sigma1278b', coverage => '45x', },
        { taxon => 643680, name => 'EC1118',     coverage => 'unknown', },

        # wustl 11 yeast strains
        { taxon => 929587, name => 'CBS_7960', coverage => '17x', },
        { taxon => 464025, name => 'CLIB215',  coverage => '16.9x', },
        { taxon => 929629, name => 'CLIB324',  coverage => '7.14x', },
        { taxon => 947035, name => 'CLIB382',  coverage => '5.96x', },
        { taxon => 947036, name => 'FL100',    coverage => '7.1x', },
        { taxon => 947039, name => 'PW5',      coverage => '16.10x', },
        { taxon => 929585, name => 'T7',       coverage => '25.4x', },
        { taxon => 471859, name => 'T73',      coverage => '13.9x', },
        { taxon => 947040, name => 'UC5',      coverage => '15.7x', },
        { taxon => 462210, name => 'Y10',      coverage => '6.6x', },
        { taxon => 929586, name => 'YJM269',   coverage => '16.7x', },

        # wine
        { taxon => 764097, name => 'AWRI796',     coverage => '20x', },
        { taxon => 764098, name => 'Lalvin_QA23', coverage => '20x', },
        { taxon => 764099, name => 'Vin13',       coverage => '20x', },
        { taxon => 764100, name => 'VL3',         coverage => '20x', },
        { taxon => 764101, name => 'FostersO',    coverage => '20x', },
        { taxon => 764102, name => 'FostersB',    coverage => '20x', },

        #{ taxon => 545124, name => 'AWRI1631',   coverage => '7x', },
        #{ taxon => 538975, name => 'M22',        coverage => '2.6x', },
        #{ taxon => 538976, name => 'YPS163',     coverage => '2.8x', },
    );

    my $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
gunzip -c [% item.name %]/*.gz > [% item.name %]/[% item.name %].fasta
perl -p -i -e '/>/ and s/\>gi\|(\d+).*/\>gi_$1/' [% item.name %]/[% item.name %].fasta
RepeatMasker [% item.name %]/*.fasta -species fungi -xsmall -s --parallel 8
mv [% item.name %]/[% item.name %].fasta.masked [% item.name %]/[% item.name %].fa
perl [% pl_dir %]/blastz/bz.pl -dt [% data_dir %]/S288C_58 -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/S288Cvs[% item.name %]_58 -s set11 -p 6 
perl [% pl_dir %]/alignDB/extra/two_way_batch.pl -d S288Cvs[% item.name %] -e yeast_58 -t="4932,S288C" -q "[% item.taxon %],[% item.name %]" -a [% data_dir %]/S288Cvs[% item.name %]_58 -at 10000 -st 1000000 --parallel 4 --run all

[% END -%]
EOF

    $tt->process( \$text,
        { data => \@data, data_dir => $data_dir, pl_dir => $pl_dir, },
        "auto_wgs.sh" )
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
cd [% data_dir %]

[% FOREACH item IN data -%]
# [% item.name %] [% item.coverage %]
RepeatMasker [% item.name %]/*.fasta -species fungi -xsmall -s --parallel 8
mv [% item.name %]/[% item.name %].fasta.masked [% item.name %]/[% item.name %].fa
perl [% pl_dir %]/blastz/bz.pl -dt [% data_dir %]/S288C_58 -dq [% data_dir %]/[% item.name %] -dl [% data_dir %]/S288Cvs[% item.name %]_58 -s set11 -p 6 
perl [% pl_dir %]/alignDB/extra/two_way_batch.pl -d S288Cvs[% item.name %] -e yeast_58 -t="4932,S288C" -q "[% item.taxon %],[% item.name %]" -a [% data_dir %]/S288Cvs[% item.name %]_58 -at 10000 -st 1000000 --parallel 4 --run all

[% END -%]
EOF

    $tt->process( \$text,
        { data => \@data, data_dir => $data_dir, pl_dir => $pl_dir, },
        "auto_sgrp.sh" )
        or die Template->error;
}

{    # multi
    my $tt   = Template->new;
    my @data = (
        {   goal_db  => 'S288CvsThree',
            dbs      => 'S288CvsSpar,S288CvsRM11,S288CvsYJM789',
            outgroup => '0query',
            target   => '0target',
            queries  => '1query,2query',
            all_freq => 3,
        },
        {   goal_db => 'S288CvsSix',
            dbs     => 'S288CvsSpar,S288CvsRM11,S288CvsYJM789'
                . ',S288CvsDBVPG6765,S288CvsSK1,S288CvsY55',
            outgroup => '0query',
            target   => '0target',
            queries  => '1query,2query,3query,4query,5query',
            all_freq => 6,
        },
        {   goal_db => 'S288CvsTen',
            dbs     => 'S288CvsSpar,S288CvsRM11,S288CvsYJM789'
                . ',S288CvsJAY291,S288CvsSigma1278b,S288CvsEC1118'
                . ',S288CvsY55,S288CvsSK1,S288CvsW303,S288CvsDBVPG6765',
            outgroup => '0query',
            target   => '0target',
            queries  => '1query,2query,3query,4query,5query'
                . ',6query,7query,8query,9query',
            all_freq => 10,
        },
    );

    my $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

[% FOREACH item IN data -%]
# [% item.goal_db %]
perl [% pl_dir %]/alignDB/extra/join_dbs.pl --dbs [% item.dbs %] --goal_db [% item.goal_db _ "_10k" %] --outgroup [% item.outgroup %] --target [% item.target %] --queries [% item.queries %] --no_insert=1 --trimmed_fasta=1 --length 10000

perl [% pl_dir %]/alignDB/extra/multi_way_batch.pl -d [% item.goal_db _ "_10k" %] -e yeast_58 -f [% data_dir %]/[% item.goal_db _ "_10k" %] --all_freq [% item.all_freq %] -lt 10000 -st 100000 --parallel=6 --run all

[% END -%]
EOF

    $tt->process( \$text,
        { data => \@data, data_dir => $data_dir, pl_dir => $pl_dir, },
        "auto_joins.sh" )
        or die Template->error;
}
