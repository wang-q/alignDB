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

    Bacillus_anthracis
    Burkholderia_mallei
    Mycobacterium_tuberculosis
    Pseudomonas_aeruginosa
    Streptococcus_agalactiae
    Vibrio_cholerae
    Xanthomonas_campestris
    Xanthomonas_oryzae
};

my @unused = qw{
    Bacillus_anthracis
    Burkholderia_mallei
    Mycobacterium_tuberculosis
    Pseudomonas_aeruginosa
    Streptococcus_agalactiae
    Vibrio_cholerae
    Xanthomonas_campestris
    Xanthomonas_oryzae
};
my %filter = map { $_ => 1 } @unused;

my @data;
for my $i ( 0 .. $#names ) {
    my $name = $names[$i];
    my $item = {
        name       => $name,
        text       => join( " ", split /_/, $name ),
        inter      => 0,
        common_file => "$name.common.chart.xlsx",
        gc_file    => "$name.gc.chart.xlsx",
    };

    push @data, $item;
}

my $tt = Template->new;
my $text;

$text = <<'EOF';
[% FOREACH item IN data -%]
[%   chart = 'common' -%]
[%  IF item.inter == 1;
        replace = '--replace diversity=divergence';
    ELSE;
        replace = '';
    END
-%]
REM [% item.name %]
if not exist [% item.name %].[% chart %].xlsx goto skip[% item.name %]
perl d:/wq/Scripts/alignDB/stat/[% chart %]_chart_factory.pl [% replace %] -i [% item.name %].[% chart %].xlsx
:skip[% item.name %]

[% END -%]
EOF
$tt->process( \$text, { data => \@data, }, "bac_common_chart.bat" )
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
perl d:/wq/Scripts/alignDB/stat/[% chart %]_chart_factory.pl [% replace %] --add_trend 1 -i [% item.name %].[% chart %].xlsx
:skip[% item.name %]

[% END -%]
EOF
$tt->process( \$text, { data => \@data, }, "bac_gc_chart.bat" )
    or die Template->error;

$text = <<'EOF';
[% USE Math -%]
---
charts:
[% FOREACH item IN data -%]
  [% item.common_file %]:
    d1_comb_pi_gc_cv:
      4:
        - [% loop.index % 4 %]
        - [% Math.int(loop.index / 4) %]
[% END -%]
texts:
[% FOREACH item IN data -%]
  - text: [% item.text %]
    size: 8
    bold: 1
    italic: 1
    pos:
      - [% loop.index % 4 %]
      - [% Math.int(loop.index / 4) - 0.05 %]
[% END -%]
EOF

#$tt->process( \$text, { data => [ grep { !$filter{ $_->{name} } } @data ], },
#    'Fig_bac_d1_gc_cv.yml' )
#    or die Template->error;

$tt->process( \$text, { data => [  @data ], },
    'Fig_bac_d1_gc_cv.yml' )
    or die Template->error;

__END__
