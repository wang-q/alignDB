#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;
use Template;
use YAML qw(Dump Load DumpFile LoadFile);

my @data = (
    {   name  => "HumanvsRhesus",
        tag   => "pair",
        text  => "Human versus Rhesus",
        inter => 1,
    },
    {   name  => "HumanvsGibbon",
        tag   => "pair",
        text  => "Human versus Gibbon",
        inter => 1,
    },
    {   name   => "MousevsXIIS",
        tag    => "multi",
        text   => "12 mouse lines",
        number => 12,
        inter  => 0,
    },
    {   name  => "MousevsSpretus_Ei",
        tag   => "pair",
        text  => "Mouse versus M. spretus",
        inter => 1,
    },
    {   name  => "DsimvsDsec",
        tag   => "pair",
        text  => "D. simulans versus D. sechellia",
        inter => 1,
    },
    {   name  => "DyakvsDere",
        tag   => "pair",
        text  => "D. yakuba versus D. erecta",
        inter => 1,
    },
    {   name   => "NipvsXXIV",
        tag    => "multi",
        text   => "24 japonica rice lines",
        number => 24,
        inter  => 0,
    },
    {   name  => "Nipvs9311",
        tag   => "pair",
        text  => "Nipponbare versus 93-11",
        inter => 0,
    },
    {   name  => "AthvsLyrata",
        tag   => "pair",
        text  => "A. thaliana versus A. lyrata",
        inter => 1,
    },
    {   name  => "S288CvsSpar",
        tag   => "pair",
        text  => "S. cerevisiae versus S. paradoxus",
        inter => 1,
    },
    {   name  => "AfumvsNfis",
        tag   => "pair",
        text  => "A. fumigatus versus N. fischeri",
        inter => 1,
    },
    {   name  => "AoryvsAfla",
        tag   => "pair",
        text  => "A. oryzae versus A. flavus",
        inter => 1,
    },

    # MousevsCAST_Ei got bad segment_cv_indel_3 result
    #{   name   => "MousevsXIIC",
    #    tag    => "multi",
    #    text   => "12 mouse lines",
    #    number => 12,
    #},
    # bad combined_distance
    #{   name  => "PbervsPcha",
    #    tag   => "pair",
    #    text  => "P. berghei versus P. chabaudi",
    #    inter => 1,
    #},
    # bad segment_cv_indel_3
    #{   name  => "MousevsRat",
    #    tag   => "pair",
    #    text  => "Mouse versus Rat",
    #    inter => 1,
    #},
);

for my $i ( 0 .. $#data ) {
    my $name = $data[$i]->{name};
    if ( $data[$i]->{tag} eq "pair" ) {
        $data[$i]->{common_file} = "$name.common.chart.xls";
        $data[$i]->{d1_sheet}    = "combined_distance";
    }
    else {
        $data[$i]->{common_file} = "$name.multi.chart.xls";
        $data[$i]->{d1_sheet}    = "combined_pigccv";
    }
    $data[$i]->{gc_file}      = "$name.gc.chart.xls";
    $data[$i]->{gene_file}    = "$name.gene.chart.xls";
    $data[$i]->{coding_sheet} = "coding_all";
}

my $tt = Template->new;
my $text;

$text = <<'EOF';
[% FOREACH item IN data -%]
[%  IF item.tag == 'pair';
        chart = 'common';
    ELSE;
        chart = 'multi';
    END
-%]
REM [% item.name %]
perl d:/wq/Scripts/alignDB/stat/[% chart %]_stat_factory.pl -s 114.212.202.159 -d [% item.name %]

[% END -%]
EOF
$tt->process( \$text, { data => \@data, }, "pair_common_stat.bat" )
    or die Template->error;

$text = <<'EOF';
[% FOREACH item IN data -%]
[%  IF item.tag == 'pair';
        chart = 'common';
    ELSE;
        chart = 'multi';
    END
-%]
[%  IF item.inter == 1;
        replace = '--replace diversity=divergence';
    ELSE;
        replace = '';
    END
-%]
REM [% item.name %]
if not exist [% item.name %].[% chart %].xlsx goto skip[% item.name %]
perl d:/wq/Scripts/alignDB/stat/[% chart %]_chart_factory.pl [% replace %] -i [% item.name %].[% chart %].xlsx
perl d:/wq/Scripts/alignDB/fig/xlsx2xls.pl -d [% item.name %].[% chart %].chart.xlsx
:skip[% item.name %]

[% END -%]
EOF
$tt->process( \$text, { data => \@data, }, "pair_common_chart.bat" )
    or die Template->error;

$text = <<'EOF';
[% FOREACH item IN data -%]
[%  chart = 'gc' -%]
[%  IF item.inter == 1;
        replace = '--replace diversity=divergence';
    ELSE;
        replace = '';
    END
-%]
REM [% item.name %]
if not exist [% item.name %].[% chart %].xlsx goto skip[% item.name %]
perl d:/wq/Scripts/alignDB/stat/[% chart %]_chart_factory.pl [% replace %] -i [% item.name %].[% chart %].xlsx
perl d:/wq/Scripts/alignDB/fig/xlsx2xls.pl -d [% item.name %].[% chart %].chart.xlsx
:skip[% item.name %]

[% END -%]
EOF
$tt->process( \$text, { data => \@data, }, "pair_gc_chart.bat" )
    or die Template->error;

$text = <<'EOF';
[% FOREACH item IN data -%]
[%  chart = 'gene' -%]
[%  IF item.inter == 1;
        replace = '--replace diversity=divergence';
    ELSE;
        replace = '';
    END
-%]
REM [% item.name %]
if not exist [% item.name %].[% chart %].xlsx goto skip[% item.name %]
perl d:/wq/Scripts/alignDB/stat/[% chart %]_chart_factory.pl [% replace %] -i [% item.name %].[% chart %].xlsx
perl d:/wq/Scripts/alignDB/fig/xlsx2xls.pl -d [% item.name %].[% chart %].chart.xlsx
:skip[% item.name %]

[% END -%]
EOF
$tt->process( \$text, { data => \@data, }, "pair_gene_chart.bat" )
    or die Template->error;

$text = <<'EOF';
[% USE Math -%]
---
charts:
[% FOREACH item IN data -%]
  [% item.common_file %]:
    [% item.d1_sheet %]:
      4:
        - [% loop.index % 2 %]
        - [% Math.int(loop.index / 2) %]
[% END -%]
texts:
[% FOREACH i IN data -%]
  - text: [% lables.shift %]
    size: 12
    bold: 1
    pos:
      - [% loop.index % 2 %]
      - [% Math.int(loop.index / 2) - 0.05 %]
[% END -%]
[% FOREACH item IN data -%]
  - text: [% item.text %]
    size: 7
    bold: 1
    pos:
      - [% loop.index % 2 + 0.1 %]
      - [% Math.int(loop.index / 2) - 0.05 %]
[% END -%]

EOF
$tt->process( \$text, { data => \@data, lables => [ 'a' .. 'z' ], },
    "Fig_pair_d1_gc_cv.yml" )
    or die Template->error;

$text = <<'EOF';
---
charts:
[% FOREACH item IN data -%]
  [% item.gc_file %]:
    segment_gc_indel_3:
      1:
        - 0
        - [% loop.index %]
      2:
        - 1
        - [% loop.index %]
    segment_cv_indel_3:
      1:
        - 2
        - [% loop.index %]
      2:
        - 3
        - [% loop.index %]
[% END -%]
texts:
[% FOREACH i IN data -%]
  - text: [% lables.shift %]
    size: 12
    bold: 1
    pos:
      - 0
      - [% loop.index - 0.05 %]
[% END -%]
[% FOREACH item IN data -%]
  - text: [% item.text %]
    size: 7
    bold: 1
    pos:
      - 0.1
      - [% loop.index - 0.05 %]
[% END -%]

EOF
$tt->process( \$text, { data => \@data, lables => [ 'a' .. 'z' ], },
    "Fig_pair_gc_cv_indel.yml" )
    or die Template->error;

__END__
