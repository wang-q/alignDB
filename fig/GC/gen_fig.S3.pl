#!/usr/bin/perl
use strict;
use warnings;

use Template;
use YAML qw(Dump Load DumpFile LoadFile);

my $text = <<'EOF';
---
charts:
[% FOREACH item IN data -%]
  [% item.mgc_file %]:
    segment_gc_indel_3:
      1:
        - 0
        - [% item.idx %]
      2:
        - 1
        - [% item.idx %]
    segment_cv_indel_3:
      1:
        - 2
        - [% item.idx %]
      2:
        - 3
        - [% item.idx %]
[% END -%]
texts:
[% FOREACH item IN data -%]
  - text: [% item.species %]
    size: 12
    bold: 1
    italic: 1
    rotate: 0
    pos:
      - 0
      - [% item.idx - 0.05 %]
[% END -%]

EOF

my $names = [
    [   qw{
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
            }
    ],
    [   qw{
            Escherichia_coli
            Francisella_tularensis
            Haemophilus_influenzae
            Helicobacter_pylori
            Lactococcus_lactis
            Legionella_pneumophila
            Listeria_monocytogenes
            Methanococcus_maripaludis
            Neisseria_meningitidis
            Prochlorococcus_marinus
            Pseudomonas_putida
            Rhodopseudomonas_palustris
            }
    ],
    [   qw{
            Salmonella_enterica
            Shewanella_baltica
            Staphylococcus_aureus
            Streptococcus_equi
            Streptococcus_pneumoniae
            Streptococcus_pyogenes
            Streptococcus_suis
            Streptococcus_thermophilus
            Sulfolobus_islandicus
            Xylella_fastidiosa
            Yersinia_pestis
            Yersinia_pseudotuberculosis
            }
    ],
];

for my $i ( 0 .. scalar @{$names} - 1 ) {
    my $output = "Fig_S3_part$i.yml";
    my @data;
    for my $j ( 0 .. scalar @{ $names->[$i] } - 1 ) {
        my $name = $names->[$i][$j];
        my $item = {
            name       => $name,
            species    => join( " ", split /_/, $name ),
            idx        => $j,
            multi_file => "$name.multi.chart.xlsx",
            mgc_file   => "$name.gc.chart.xlsx",
        };

        push @data, $item;
    }

    my $tt = Template->new;
    $tt->process( \$text, { data => \@data, }, $output ) or die Template->error;
}

__END__

