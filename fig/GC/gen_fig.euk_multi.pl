#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;
use Template;
use YAML qw(Dump Load DumpFile LoadFile);

my @data = (
    {   name   => "HumanvsCGOR",
        tag    => "multi",
        text   => "Human, chimp, gorilla and orangutan",
        number => 4,
        inter  => 1,
    },

    # MousevsCAST_Ei got bad segment_cv_indel_3 result
    #{   name   => "MousevsXIIC",
    #    tag    => "multi",
    #    text   => "12 mouse lines",
    #    number => 12,
    #},
    {   name   => "MousevsXIIS",
        tag    => "multi",
        text   => "12 mouse lines",
        number => 12,
        inter  => 0,
    },
    {   name   => "DmelvsXXII",
        tag    => "multi",
        text   => "22 Drosophila lines",
        number => 22,
        inter  => 0,
    },
    {   name   => "NipvsXXIV",
        tag    => "multi",
        text   => "24 japonica rice lines",
        number => 24,
        inter  => 0,
    },
    {   name   => "AthvsXIX",
        tag    => "multi",
        text   => "19 Arabidopsis lines",
        number => 19,
        inter  => 0,
    },
    {   name   => "S288CvsVIII_WGS",
        tag    => "multi",
        text   => "8 yeast lines",
        number => 8,
        inter  => 0,
    },
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
$tt->process( \$text, { data => \@data, }, "multi_common_stat.bat" )
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
$tt->process( \$text, { data => \@data, }, "multi_common_chart.bat" )
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
$tt->process( \$text, { data => \@data, }, "multi_gc_chart.bat" )
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
$tt->process( \$text, { data => \@data, }, "multi_gene_chart.bat" )
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
    "Fig_multi_d1_gc_cv.yml" )
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
    "Fig_multi_gc_cv_indel.yml" )
    or die Template->error;

$text = <<'EOF';
---
charts:
[% FOREACH item IN data -%]
  [% item.common_file %]:
    pigccv_freq_1of[% item.number %]:
      4:
        - 0
        - [% loop.index %]
    pigccv_freq_2of[% item.number %]:
      4:
        - 1
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
  - text: 1/[% item.number %]
    size: 7
    bold: 1
    pos:
      - 0.4
      - [% loop.index + 0.06 %]
  - text: 2/[% item.number %]
    size: 7
    bold: 1
    pos:
      - 1.4
      - [% loop.index + 0.06 %]
[% END -%]

EOF
$tt->process( \$text, { data => \@data, lables => [ 'a' .. 'z' ], },
    "Fig_gc_freq.yml" )
    or die Template->error;

$text = <<'EOF';
---
charts:
[% FOREACH item IN data -%]
  [% item.gc_file %]:
    combined_distance:
      2:
        - 0
        - [% loop.index %]
    combined_a2d:
      2:
        - 1
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
    "Fig_gc_trough.yml" )
    or die Template->error;

$text = <<'EOF';
---
charts:
[% FOREACH item IN data -%]
  [% item.gene_file %]:
    [% item.coding_sheet %]:
      1:
        - 0
        - [% loop.index %]
      2:
        - 1
        - [% loop.index %]
      3:
        - 2
        - [% loop.index %]
      4:
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
    "Fig_coding_border.yml" )
    or die Template->error;

__END__
