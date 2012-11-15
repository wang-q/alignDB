#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use Template;

my @dbs = qw{
    HumanvsCGOR
    MousevsXIIS
    DmelvsXXII
    NipvsXXIV
    AthvsXIX
    S288CvsVIII_WGS
    
    HumanvsGibbon
    HumanvsRhesus
    DsimvsDsec
    DyakvsDere
    Nipvs9311
    AthvsLyrata
    S288CvsSpar
    AfumvsNfis
    AoryvsAfla
    PbervsPcha

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
    Methanococcus_maripaludis
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
    Sulfolobus_islandicus
    Xylella_fastidiosa
    Yersinia_pestis
    Yersinia_pseudotuberculosis
};

my $tt = Template->new;
my $text;

{
    $text = <<'EOF';
#!/bin/bash

[% FOREACH db IN dbs -%]
perl /home/wangq/Scripts/alignDB/util/dup_db.pl -d [% db %] -f /home/wangq/data/GC/[% db %].sql
[% END -%]

EOF

    $tt->process( \$text, { dbs => \@dbs, }, "gc_dump.sh" )
        or die Template->error;
}

{
    $text = <<'EOF';
#!/bin/bash

[% FOREACH db IN dbs -%]
perl /home/wangq/Scripts/alignDB/util/dup_db.pl -f /home/wangq/data/GC/[% db %].sql -g [% db %] 
[% END -%]

EOF

    $tt->process( \$text, { dbs => \@dbs, }, "gc_restore.sh" )
        or die Template->error;
}
