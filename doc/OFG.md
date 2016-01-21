# OFG

## spo11 (cell 2011)

### Get data

Original paper:

Pan, J. et al. A Hierarchical Combination of Factors Shapes the Genome-wide Topography of Yeast Meiotic Recombination Initiation. Cell 144, 719Ð731 (2011).

```bash
mkdir -p ~/data/ofg/spo11
cd ~/data/ofg/spo11

wget -N http://www.cell.com/cms/attachment/612087/4902311/mmc2.xls

perl ~/Scripts/fig_table/xlsx2csv.pl -f mmc2.xls --sheet 'Spo11 Hot Spot Annotations' \
    | perl -F/,/ -an -e '
    /^CHROM/ and next;
    /^chr/ or next;
    $F[0] =~ s/^chr//;
    print qq{$F[0]\t$F[1]\t$F[2]\n};
    ' \
    > spo11_hot.bed
```

### Process

1. S288Cvsself center_intact

    ```bash
    cd ~/data/ofg/spo11

    perl ~/Scripts/alignDB/util/dup_db.pl -g S288cvsself_spo11 -f ~/data/dumps/mysql/S288cvsself.sql.gz

    perl ~/Scripts/alignDB/ofg/insert_bed.pl \
        -d S288cvsself_spo11 \
        --style center_intact --batch 1 --parallel 8 \
        --dG \
        --tag spo11 --type hotspot -f ~/data/ofg/spo11/spo11_hot.bed

    perl ~/Scripts/alignDB/init/update_sw_cv.pl -d S288cvsself_spo11 --batch 1 --parallel 8
    perl ~/Scripts/alignDB/init/update_feature.pl \
        -d S288cvsself_spo11 \
        -e saccharomyces_cerevisiae_core_29_82_4 \
        --batch 1 --parallel 8

    perl ~/Scripts/alignDB/stat/ofg_stat_factory.pl -d S288cvsself_spo11 --index --chart
    ```

2. Scer_n8_pop center_intact

    ```bash
    cd ~/data/ofg/spo11

    perl ~/Scripts/alignDB/util/dup_db.pl -g Scer_n8_pop_spo11 -f ~/data/dumps/mysql/Scer_n8_pop.sql.gz

    perl ~/Scripts/alignDB/ofg/insert_bed.pl \
        -d Scer_n8_pop_spo11 \
        --style center_intact --batch 10 --parallel 8 \
        --tag spo11 --type hotspot -f ~/data/ofg/spo11/spo11_hot.bed

    perl ~/Scripts/alignDB/init/update_sw_cv.pl -d Scer_n8_pop_spo11 --batch 10 --parallel 8
    perl ~/Scripts/alignDB/init/update_feature.pl \
        -d Scer_n8_pop_spo11 \
        -e saccharomyces_cerevisiae_core_29_82_4 \
        --batch 10 --parallel 8

    perl ~/Scripts/alignDB/stat/ofg_stat_factory.pl -d Scer_n8_pop_spo11 --index --chart
    ```

3. S288Cvsself edge

    ```bash
    cd ~/data/ofg/spo11

    perl ~/Scripts/alignDB/util/dup_db.pl -g S288cvsself_spo11_edge -f ~/data/dumps/mysql/S288cvsself.sql.gz

    perl ~/Scripts/alignDB/ofg/insert_bed.pl \
        -d S288cvsself_spo11_edge \
        --style edge --batch 1 --parallel 8 \
        --tag spo11 --type hotspot -f ~/data/ofg/spo11/spo11_hot.bed

    perl ~/Scripts/alignDB/init/update_sw_cv.pl -d S288cvsself_spo11_edge --batch 1 --parallel 8
    perl ~/Scripts/alignDB/init/update_feature.pl \
        -d S288cvsself_spo11_edge \
        -e saccharomyces_cerevisiae_core_29_82_4 \
        --batch 1 --parallel 8

    perl ~/Scripts/alignDB/stat/ofg_stat_factory.pl -d S288cvsself_spo11_edge --index --chart
    ```
