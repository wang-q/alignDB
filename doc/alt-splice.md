# alt-splice

## Flanking introns

### Put data to `~/data/alt-splice/ase_flanking/`.

```
$ perl ~/Scripts/tool/list_dirtree.pl ~/data/alt-splice/ase_flanking/
+-[ase_flanking]
    |-ASase.bed                                   |324240 lines |     6.46M
    |-AScse.bed                                   | 78702 lines |     1.57M
    |-nonAS.bed                                   |  9296 lines |    189.8K
```

### Processing

```bash
cd ~/data/alt-splice/ase_flanking/

# Runtime 29 minutes and 53 seconds.
perl ~/Scripts/alignDB/util/dup_db.pl -f ~/data/dumps/mysql/Human_n11cg_chimp_basic.sql.gz -g Human_n11cg_chimp_as_flanking

perl ~/Scripts/alignDB/ofg/insert_bed.pl \
    -d Human_n11cg_chimp_as_flanking \
    --style center_intact \
    --parallel 8 \
    --tag intron --type ASase -f ~/data/alt-splice/ase_flanking/ASase.bed \
    --tag intron --type AScse -f ~/data/alt-splice/ase_flanking/AScse.bed \
    --tag intron --type nonAS -f ~/data/alt-splice/ase_flanking/nonAS.bed

perl ~/Scripts/alignDB/init/update_sw_cv.pl \
    -d Human_n11cg_chimp_as_flanking \
    --parallel 8

perl ~/Scripts/alignDB/init/update_feature.pl \
    -d Human_n11cg_chimp_as_flanking \
    -e homo_sapiens_core_82_37 \
    --parallel 8

perl /home/wangq/Scripts/alignDB/stat/ofg_stat_factory.pl \
    --by type -d Human_n11cg_chimp_as_flanking \
    -o Human_n11cg_chimp_as_flanking.ofg.xlsx

```
