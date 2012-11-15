#!/usr/bin/perl
use strict;
use warnings;

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
    
    # bad combined_distance
    #{   name  => "PbervsPcha",
    #    tag   => "pair",
    #    text  => "P. berghei versus P. chabaudi",
    #    inter => 1,
    #},
    
    {   name  => "MousevsSpretus_Ei",
        tag   => "pair",
        text  => "Mouse versus M. spretus",
        inter => 1,
    },
    
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
        $data[$i]->{common_file} = "$name.common.chart.xlsx";
    }
    else {
        $data[$i]->{common_file} = "$name.multi.chart.xlsx";
    }
    $data[$i]->{gc_file}      = "$name.gc.chart.xlsx";
}

my $tt = Template->new;
my $text;

$text = <<'EOF';
---
autofit: A:J
texts:
  - text: "Species"
    pos: A2:A3
    merge: 1
  - text: "No. of lines"
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
  - text: "SNPs Per 100 bp"
    pos: F2:F3
    merge: 1
  - text: "Correlation coefficients (r) between"
    pos: G2:J2
    merge: 1
  - text: "GC & Pi"
    pos: G3
  - text: "GC & Indel"
    pos: H3
  - text: "CV & Pi"
    pos: I3
  - text: "CV & Indel"
    pos: J3
  - text: "r_squared"
    pos: K2:N2
    merge: 1
  - text: "GC & Pi"
    pos: K3
  - text: "GC & Indel"
    pos: L3
  - text: "CV & Pi"
    pos: M3
  - text: "CV & Indel"
    pos: N3
  - text: "p_value"
    pos: O2:R2
    merge: 1
  - text: "GC & Pi"
    pos: O3
  - text: "GC & Indel"
    pos: P3
  - text: "CV & Pi"
    pos: Q3
  - text: "CV & Indel"
    pos: R3
  - text: "slope"
    pos: S2:V2
    merge: 1
  - text: "GC & Pi"
    pos: S3
  - text: "GC & Indel"
    pos: T3
  - text: "CV & Pi"
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
  [% item.common_file %]:
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
  [% item.gc_file %]:
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
$tt->process( \$text, { data => \@data, }, 'Table_S_euk_pair.yml' )
    or die Template->error;

__END__
