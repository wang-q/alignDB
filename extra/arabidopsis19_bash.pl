#!/usr/bin/perl
use strict;
use warnings;

use Template;
use File::Basename;
use File::Find::Rule;
use File::Remove qw(remove);
use File::Spec;
use String::Compare;
use YAML qw(Dump Load DumpFile LoadFile);

my $store_dir = shift || "/home/wangq/data/alignment/arabidopsis19";
{    # on linux
    my $data_dir    = "/home/wangq/data/alignment/arabidopsis19";
    my $pl_dir      = "/home/wangq/Scripts";
    my $kentbin_dir = "/home/wangq/bin/x86_64";

    # nature 2011
    my $seq_dir = "/home/wangq/data/1001/19genomes/fasta/MASKED";

    my $tt = Template->new;

    my @data = (
        { taxon => 900201, name => "Bur_0",  origin => "Ireland" },
        { taxon => 900202, name => "Can_0",  origin => "Canary Isles" },
        { taxon => 900203, name => "Ct_1",   origin => "Italy" },
        { taxon => 900204, name => "Edi_0",  origin => "Scotland" },
        { taxon => 900205, name => "Hi_0",   origin => "Netherlands" },
        { taxon => 900206, name => "Kn_0",   origin => "Lithuania" },
        { taxon => 900207, name => "Ler_0",  origin => "Poland" },
        { taxon => 900208, name => "Mt_0",   origin => "Libya" },
        { taxon => 900209, name => "No_0",   origin => "Germany" },
        { taxon => 900210, name => "Oy_0",   origin => "Norway" },
        { taxon => 900211, name => "Po_0",   origin => "Germany" },
        { taxon => 900212, name => "Rsch_4", origin => "Russia" },
        { taxon => 900213, name => "Sf_2",   origin => "Spain" },
        { taxon => 900214, name => "Tsu_0",  origin => "Japan" },
        { taxon => 900215, name => "Wil_2",  origin => "Russia" },
        { taxon => 900216, name => "Ws_0",   origin => "Russia" },
        { taxon => 900217, name => "Wu_0",   origin => "Germany" },
        { taxon => 900218, name => "Zu_0",   origin => "Germany" },
    );

    my @files = File::Find::Rule->file->name('*.fas')->in($seq_dir);

    for my $item ( sort @data ) {

        # match the most similar name
        my ($file) = map { $_->[0] }
            sort { $b->[1] <=> $a->[1] }
            map { [ $_, compare( lc basename($_), lc $item->{name} ) ] } @files;
        $item->{file} = $file;

        # prepare working dir
        my $dir = File::Spec->catdir( $data_dir, $item->{name} );
        mkdir $dir if !-e $dir;
        $item->{dir} = $dir;
    }

    my $basecount = File::Spec->catfile( $data_dir, "basecount.txt" );
    remove( \1, $basecount ) if -e $basecount;

    # taxon.csv
    my $text = <<'EOF';
[% FOREACH item IN data -%]
[% item.taxon %],Arabidopsis,thaliana,[% item.name %],,
[% END -%]
EOF
    $tt->process(
        \$text,
        { data => \@data, },
        File::Spec->catfile( $store_dir, "taxon.csv" )
    ) or die Template->error;

    # chr_length.csv
    $text = <<'EOF';
[% FOREACH item IN data -%]
[% item.taxon %],chrUn,999999999,[% item.name %]/arabidopsis19/1001
[% END -%]
EOF
    $tt->process(
        \$text,
        { data => \@data, },
        File::Spec->catfile( $store_dir, "chr_length.csv" )
    ) or die Template->error;

    #
    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# basecount and split
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.origin %]
echo [% item.name %] >> [% data_dir %]/basecount.txt
[% kentbin_dir %]/faCount [% item.file %] >> [% data_dir %]/basecount.txt
echo >> [% data_dir %]/basecount.txt

[% kentbin_dir %]/faSplit byname [% item.file %] [% item.dir %]/

find [% item.dir %] -name "*.fa" | sed "s/\.fa$//" | xargs -i echo mv {}.fa {}.fasta | sh

[% END -%]

#----------------------------#
# repeatmasker
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.origin %]
RepeatMasker [% item.dir %]/*.fasta -species arabidopsis -xsmall -s --parallel 4

find [% item.dir %] -name "*.fasta.masked" | sed "s/\.fasta\.masked$//" | xargs -i echo mv {}.fasta.masked {}.fa | sh

[% END -%]

EOF

    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "auto_ath19_file.sh" )
    ) or die Template->error;

    $text = <<'EOF';
#!/bin/bash
cd [% data_dir %]

#----------------------------#
# blastz
#----------------------------#
[% FOREACH item IN data -%]
# [% item.name %] [% item.origin %]
perl [% pl_dir %]/blastz/bz.pl -dt [% data_dir %]/ath_65 -dq [% data_dir %]/[% item.name %] \
    -dl [% data_dir %]/Athvs[% item.name %] -s set01 -p 4 --noaxt -pb lastz --lastz --paired

perl [% pl_dir %]/blastz/lpcna.pl -dt [% data_dir %]/ath_65 -dq [% data_dir %]/[% item.name %] \
    -dl [% data_dir %]/Athvs[% item.name %]

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data        => \@data,
            data_dir    => $data_dir,
            pl_dir      => $pl_dir,
            kentbin_dir => $kentbin_dir
        },
        File::Spec->catfile( $store_dir, "auto_ath19_bz.sh" )
    ) or die Template->error;
}

{    # on windows
    my $data_dir = "d:/data/alignment/arabidopsis19";
    my $pl_dir   = "d:/wq/Scripts";

    my $tt = Template->new;

    my @data = (
        { taxon => 900201, name => "Bur_0",  coverage => 25, },
        { taxon => 900202, name => "Can_0",  coverage => 47, },
        { taxon => 900203, name => "Ct_1",   coverage => 50, },
        { taxon => 900204, name => "Edi_0",  coverage => 52, },
        { taxon => 900205, name => "Hi_0",   coverage => 33, },
        { taxon => 900206, name => "Kn_0",   coverage => 28, },
        { taxon => 900207, name => "Ler_0",  coverage => 27, },
        { taxon => 900208, name => "Mt_0",   coverage => 30, },
        { taxon => 900209, name => "No_0",   coverage => 38, },
        { taxon => 900210, name => "Oy_0",   coverage => 54, },
        { taxon => 900211, name => "Po_0",   coverage => 41, },
        { taxon => 900212, name => "Rsch_4", coverage => 38, },
        { taxon => 900213, name => "Sf_2",   coverage => 40, },
        { taxon => 900214, name => "Tsu_0",  coverage => 48, },
        { taxon => 900215, name => "Wil_2",  coverage => 40, },
        { taxon => 900216, name => "Ws_0",   coverage => 33, },
        { taxon => 900217, name => "Wu_0",   coverage => 26, },
        { taxon => 900218, name => "Zu_0",   coverage => 31, },
    );

    my $text = <<'EOF';
cd /d [% data_dir %]

REM #----------------------------#
REM # stat
REM #----------------------------#
[% FOREACH item IN data -%]
REM # [% item.name %] [% item.coverage %]
perl [% pl_dir %]\alignDB\extra\two_way_batch.pl -d Athvs[% item.name %] -t="3702,Ath" -q "[% item.taxon %],[% item.name %]" -a [% data_dir %]\Athvs[% item.name %] -at 10000 -st 1000000 --parallel 4 --run 1-3,21,40

[% END -%]

EOF
    $tt->process(
        \$text,
        {   data     => \@data,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
        },
        File::Spec->catfile( $store_dir, "auto_ath19_stat.bat" )
    ) or die Template->error;

    $text = <<'EOF';
cd /d [% data_dir %]

REM #----------------------------#
REM # multi
REM #----------------------------#
perl [% pl_dir %]/alignDB/extra/join_dbs.pl --dbs [% dbs %] --goal_db [% goal_db %] --outgroup [% outgroup %] --target [% target %] --queries [% queries %] --no_insert=1 --trimmed_fasta=1 --length 1000

perl [% pl_dir %]/alignDB/extra/multi_way_batch.pl -d [% goal_db %] -e ath_65 -f [% data_dir %]/[% goal_db %]  -lt 1000 -st 100000 --parallel 4 --run all

EOF

    my @names = ( "Lyrata", map { $_->{name} } @data );
    my $dbs = join ',', map { "Athvs" . $_ } @names;
    my $queries = join ',', map { $_ . "query" } ( 1 .. scalar @names - 1 );
    $tt->process(
        \$text,
        {   goal_db  => "AthvsNineteen",
            outgroup => '0query',
            target   => '0target',
            dbs      => $dbs,
            queries  => $queries,
            data_dir => $data_dir,
            pl_dir   => $pl_dir,
        },
        File::Spec->catfile( $store_dir, "auto_ath19_multi.bat" )
    ) or die Template->error;
}

#{    # multi
#    my $tt         = Template->new;
#    my $strains_of = {
#        S288CvsYJM789refSpar => [qw{ Spar YJM789 }],
#        S288CvsThree         => [qw{ Spar RM11 YJM789 }],
#        S288CvsSix           => [qw{ Spar RM11 YJM789 DBVPG6765 SK1 Y55 }],
#        S288CvsGE10M18       => [
#            qw{ Spar RM11 YJM789 JAY291 Sigma1278b EC1118 T7 AWRI796
#                Lalvin_QA23 Vin13 VL3 FostersO FostersB Kyokai_no__7 DBVPG6765
#                SK1 Y55 W303
#                }
#        ],
#        S288CvsALL32 => [
#            qw{ Spar RM11 YJM789 JAY291 Sigma1278b EC1118 CBS_7960 CLIB215
#                CLIB324 FL100 Y10 YJM269 CLIB382 PW5 T7 T73 UC5 AWRI796
#                Lalvin_QA23 Vin13 VL3 FostersO FostersB EC9_8 Kyokai_no__7
#                AWRI1631 M22 YPS163 DBVPG6765 SK1 Y55 W303
#                }
#        ],
#    };
#
#    my @data;
#    for my $dbname ( sort keys %{$strains_of} ) {
#        my @strains = @{ $strains_of->{$dbname} };
#        my $dbs     = join ',', map {"S288Cvs$_"} @strains;
#        my $queries = join ',',
#            map { $_ . "query" } ( 1 .. scalar @strains - 1 );
#        push @data,
#            {
#            goal_db  => $dbname,
#            outgroup => '0query',
#            target   => '0target',
#            dbs      => $dbs,
#            queries  => $queries,
#            };
#    }
#
#    my $text = <<'EOF';
##!/bin/bash
#cd [% data_dir %]
#
#[% FOREACH item IN data -%]
## [% item.goal_db %]
#perl [% pl_dir %]/alignDB/extra/join_dbs.pl --dbs [% item.dbs %] --goal_db [% item.goal_db %] --outgroup [% item.outgroup %] --target [% item.target %] --queries [% item.queries %] --no_insert=1 --trimmed_fasta=1 --length 1000
#
#perl [% pl_dir %]/alignDB/extra/multi_way_batch.pl -d [% item.goal_db %] -e yeast_58 -f [% data_dir %]/[% item.goal_db %]  -lt 1000 -st 100000 --parallel 4 --run all
#
#[% END -%]
#EOF
#
#    $tt->process( \$text,
#        { data => \@data, data_dir => $data_dir, pl_dir => $pl_dir, },
#        "auto_joins.sh" )
#        or die Template->error;
#}
