#!/usr/bin/perl
use strict;
use warnings;

use Template;
use YAML qw(Dump Load DumpFile LoadFile);

my @names = qw{
    Acinetobacter_baumannii
    Actinobacillus_pleuropneumoniae
    Bacillus_cereus
    Bifidobacterium_longum
    Burkholderia_cenocepacia
    Burkholderia_pseudomallei
    Campylobacter_jejuni
    Chlamydia_trachomatis
    Clostridium_botulinum
    Coxiella_burnetii
    Escherichia_coli
    Francisella_tularensis
    Haemophilus_influenzae
    Helicobacter_pylori
    Lactococcus_lactis
    Legionella_pneumophila
    Listeria_monocytogenes
    Neisseria_meningitidis
    Prochlorococcus_marinus
    Pseudomonas_putida
    Rhodopseudomonas_palustris
    Salmonella_enterica
    Shewanella_baltica
    Staphylococcus_aureus
    Streptococcus_equi
    Streptococcus_pneumoniae
    Streptococcus_pyogenes
    Streptococcus_suis
    Streptococcus_thermophilus
    Xylella_fastidiosa
    Yersinia_pestis
    Yersinia_pseudotuberculosis
    Methanococcus_maripaludis
    Sulfolobus_islandicus

    FILTERED
    Bacillus_anthracis
    Burkholderia_mallei
    Mycobacterium_tuberculosis
    Pseudomonas_aeruginosa
    Streptococcus_agalactiae
    Vibrio_cholerae
    Xanthomonas_campestris
    Xanthomonas_oryzae
};

my @data;
for my $i ( 0 .. $#names ) {
    my $name = $names[$i];
    my $item = {
        name       => $name,
        text       => join( " ", split /_/, $name ),
        tag        => 'multi',
        multi_file => "$name.multi.chart.xlsx",
        gc_file    => "$name.gc.chart.xlsx",
    };

    push @data, $item;
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
  - text: "No. of Strains"
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
  [% item.multi_file %]:
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
$tt->process( \$text, { data => \@data, }, 'Table_S_bac.yml' )
    or die Template->error;

__END__
