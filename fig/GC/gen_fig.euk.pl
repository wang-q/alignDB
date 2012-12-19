#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;
use Template;
use YAML qw(Dump Load DumpFile LoadFile);

# data in manuscript
my @data_main = (
    {   name   => "HumanvsCGOR",
        tag    => "multi",
        text   => "Human, chimp, gorilla and orangutan",
        number => 4,
        inter  => 1,
    },
    {   name   => "DmelvsXXII",
        tag    => "multi",
        text   => "22 Drosophila lines",
        number => 22,
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

# data in supplement
my @data_supp = (
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

my @data = ( @data_main, @data_supp );

for my $i ( 0 .. $#data ) {
    my $name = $data[$i]->{name};
    if ( $data[$i]->{tag} eq "pair" ) {
        $data[$i]->{common_file} = "$name.common.chart.xls";
        $data[$i]->{d1_sheet}    = "combined_distance";
        $data[$i]->{table_file} = "$name.common.chart.xlsx";
    }
    else {
        $data[$i]->{common_file} = "$name.multi.chart.xls";
        $data[$i]->{d1_sheet}    = "combined_pigccv";
        $data[$i]->{table_file} = "$name.multi.chart.xlsx";
    }
    $data[$i]->{gc_file}      = "$name.gc.chart.xls";
    $data[$i]->{gene_file}    = "$name.gene.chart.xls";
    $data[$i]->{coding_sheet} = "coding_all";
    $data[$i]->{table_gc_file}      = "$name.gc.chart.xlsx";
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
$tt->process( \$text, { data => \@data, }, "euk_common_stat.bat" )
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
$tt->process( \$text, { data => \@data, }, "euk_common_chart.bat" )
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
$tt->process( \$text, { data => \@data, }, "euk_gc_chart.bat" )
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
$tt->process( \$text, { data => \@data, }, "euk_gene_chart.bat" )
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
$tt->process( \$text, { data => \@data_main, lables => [ 'a' .. 'z' ], },
    "Fig_main_d1_gc_cv.yml" )
    or die Template->error;
$tt->process( \$text, { data => \@data_supp, lables => [ 'a' .. 'z' ], },
    "Fig_supp_d1_gc_cv.yml" )
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
$tt->process( \$text, { data => \@data_main, lables => [ 'a' .. 'z' ], },
    "Fig_main_gc_cv_indel.yml" )
    or die Template->error;
$tt->process( \$text, { data => \@data_supp, lables => [ 'a' .. 'z' ], },
    "Fig_supp_gc_cv_indel.yml" )
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
$tt->process( \$text, { data => \@data_main, lables => [ 'a' .. 'z' ], },
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
$tt->process( \$text, { data => \@data_main, lables => [ 'a' .. 'z' ], },
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
$tt->process( \$text, { data => \@data_main, lables => [ 'a' .. 'z' ], },
    "Fig_coding_border.yml" )
    or die Template->error;

$text = <<'EOF';
---
autofit: A:J
texts:
  - text: "Species"
    pos: A2:A3
    merge: 1
  - text: "No. of Genomes"
    pos: B2:B3
    merge: 1
  - text: "GC-content"
    pos: C2:C3
    merge: 1
  - text: "Aligned length (Mb)"
    pos: D2:D3
    merge: 1
  - text: "Indels Per 100 bp"
    pos: E2:E3
    merge: 1
  - text: "Substitutions Per 100 bp"
    pos: F2:F3
    merge: 1
  - text: "Correlation coefficients (r) between"
    pos: G2:J2
    merge: 1
  - text: "GC & D"
    pos: G3
  - text: "GC & Indel"
    pos: H3
  - text: "CV & D"
    pos: I3
  - text: "CV & Indel"
    pos: J3
  - text: "r_squared"
    pos: K2:N2
    merge: 1
  - text: "GC & D"
    pos: K3
  - text: "GC & Indel"
    pos: L3
  - text: "CV & D"
    pos: M3
  - text: "CV & Indel"
    pos: N3
  - text: "p_value"
    pos: O2:R2
    merge: 1
  - text: "GC & D"
    pos: O3
  - text: "GC & Indel"
    pos: P3
  - text: "CV & D"
    pos: Q3
  - text: "CV & Indel"
    pos: R3
  - text: "slope"
    pos: S2:V2
    merge: 1
  - text: "GC & D"
    pos: S3
  - text: "GC & Indel"
    pos: T3
  - text: "CV & D"
    pos: U3
  - text: "CV & Indel"
    pos: V3
  - text: =CONCATENATE(IF(S4<0,"-",""),ROUND(SQRT(K4),3),IF(O4<0.01,"**",IF(O4<0.05,"*","")))
    pos: G4
[% FOREACH item IN data -%]
  - text: [% item.text %]
    pos: A[% loop.index + 4 %]
[% END -%]
borders:
  - range: A2:J2
    top: 1
  - range: A3:J3
    bottom: 1
  - range: G2:J2
    bottom: 1
ranges:
[% FOREACH item IN data -%]
  [% item.table_file %]:
    basic:
      - copy: B2
        paste: B[% loop.index + 4 %]
      - copy: B8
        paste: C[% loop.index + 4 %]
      - copy: B4
        paste: D[% loop.index + 4 %]
      - copy: B5
        paste: E[% loop.index + 4 %]
      - copy: B6
        paste: F[% loop.index + 4 %]
  [% item.table_gc_file %]:
    segment_gc_indel_3:
      - copy: Q3
        paste: K[% loop.index + 4 %]
      - copy: Q4
        paste: O[% loop.index + 4 %]
      - copy: Q6
        paste: S[% loop.index + 4 %]
      - copy: Q18
        paste: L[% loop.index + 4 %]
      - copy: Q19
        paste: P[% loop.index + 4 %]
      - copy: Q21
        paste: T[% loop.index + 4 %]
    segment_cv_indel_3:
      - copy: Q3
        paste: M[% loop.index + 4 %]
      - copy: Q4
        paste: Q[% loop.index + 4 %]
      - copy: Q6
        paste: U[% loop.index + 4 %]
      - copy: Q18
        paste: N[% loop.index + 4 %]
      - copy: Q19
        paste: R[% loop.index + 4 %]
      - copy: Q21
        paste: V[% loop.index + 4 %]
[% END -%]
EOF
$tt->process( \$text, { data => \@data, }, 'Table_S_euk.yml' )
    or die Template->error;

__END__
