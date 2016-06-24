# Build alignments of several strains of *Saccharomyces cerevisiae*

As example for egaz and alignDB. Extract from Scer_wgs of [`OPs-download.md`](https://github.com/wang-q/withncbi/blob/master/pop/OPs-download.md) and [`OPs-align.md`](https://github.com/wang-q/withncbi/blob/master/pop/OPs-align.md).

## Build local ensembl database

Use `build_ensembl.pl`.

If don't use scripts in `gene/*`, local ensembl database names are arbitrary.

```bash
perl ~/Scripts/withncbi/ensembl/build_ensembl.pl --initdb --db saccharomyces_cerevisiae_core_29_82_4 \
    --ensembl ~/data/ensembl82/mysql/saccharomyces_cerevisiae_core_29_82_4
    
perl ~/Scripts/withncbi/ensembl/build_ensembl.pl --initdb --db yeast \
    --ensembl ~/data/ensembl82/mysql/saccharomyces_cerevisiae_core_29_82_4

```

## Two-way batch

`-chr` is omitted because all chromosomes existed in default one.

`--ensembl yeast` means different in step 20 and 22. In step 20, `yeast` is the mysql database name.
And in step 22, `yeast` is an alias to `saccharomyces_cerevisiae_core_29_82_4`.

```bash
# S288cvsRM11_1a
perl ~/Scripts/alignDB/util/two_way_batch.pl \
    -t S288c -q RM11_1a \
    -d S288cvsRM11_1a \
    -da ~/data/alignment/example/scer/Pairwise/S288cvsRM11_1a \
    -e yeast \
    -lt 5000 \
    --parallel 8 \
    -r all

```

## Multi-way batch

```bash
perl ~/Scripts/alignDB/util/multi_way_batch.pl \
    -d ScervsRM11_1a_Spar \
    -da ~/data/alignment/example/scer/Scer_n2_Spar_refined \
    --ensembl yeast \
    --outgroup \
    -lt 1000 --parallel 8 --batch 5 \
    --run all
```

## Slicing

### Yeast intergenic regions

* multi-way (fas)

```bash
mkdir -p ~/data/alignment/example/feature
cd ~/data/alignment/example/feature
 
perl ~/Scripts/alignDB/util/write_runlist_feature.pl \
    -e yeast --feature intergenic

find ~/data/alignment/example/scer/Scer_n2_Spar_refined -name "*.fas" -or -name "*.fas.gz" \
    | parallel --no-run-if-empty -j 8 \
        fasops slice {} yeast.intergenic.yml \
        --name S288c \
        -o ~/data/alignment/example/feature/{/}.fas

perl ~/Scripts/alignDB/util/multi_way_batch.pl \
    -d S288cvsRM11_1a_intergenic \
    -da . \
    -lt 1000 \
    --parallel 8 \
    --run basic
```

## Speed test

### Two-way

```bash
perl ~/Scripts/alignDB/util/two_way_batch.pl \
    -d S288cvsRM11_1a \
    -t S288c -q RM11_1a \
    -a ~/data/alignment/Ensembl/S288c/anno.yml \
    -da ~/data/alignment/example/scer/Pairwise/S288cvsRM11_1a \
    --chr ~/data/alignment/example/scer/chr_length.csv \
    -e saccharomyces_cerevisiae_core_29_82_4 \
    -lt 10000 \
    --parallel 4 --batch 10 \
    --run common
```

* On Hackintosh (i7-6700k 32G SSD)
    * `--parallel 4`: 1m44s
    * `--parallel 8`: 1m39s

* On server (E5-2690 v3 128G SSD)
    * `--parallel 4`: 2m8s
    * `--parallel 8`: 1m28s
    * `--parallel 16`: 1m15s

* On Hackintosh (i7-4790k 16G SSD)
    * `--parallel 4`: 2m52s
    * `--parallel 8`: 2m50s

* On server (E5-2690 256G HDD)
    * `--parallel 4`: 4m32s

* On macbook pro (mid 2014)
    * `--parallel 4`: 4m5s

* On desktop pc (E3-1245 v2 Windows)
    * `--parallel 4`: 6m9s

### Multi-way

```bash
perl ~/Scripts/alignDB/util/multi_way_batch.pl \
    -d Scer_n4 \
    -a ~/data/alignment/Ensembl/S288c/anno.yml \
    -da ~/data/alignment/example/scer/plan_ALL_refined \
    --chr ~/data/alignment/example/scer/chr_length.csv \
    -e saccharomyces_cerevisiae_core_29_82_4 \
    -lt 5000 \
    --parallel 4 --batch 10 \
    --run common
```

* On Hackintosh (i7-6700k 32G SSD)
    * `--parallel 4`: 5m57s
    * `--parallel 8`: 6m20s

* On server (E5-2690 v3 128G SSD)
    * `--parallel 4`: 6m2s
    * `--parallel 8`: 5m33s

### Multi-way (outgroup)

```bash
perl ~/Scripts/alignDB/util/multi_way_batch.pl \
    -d Scer_n3_Spar \
    -a ~/data/alignment/Ensembl/S288c/anno.yml \
    -da ~/data/alignment/example/scer/Scer_n3_Spar_refined \
    --chr ~/data/alignment/example/scer/chr_length.csv \
    -e saccharomyces_cerevisiae_core_29_82_4 \
    -lt 5000 \
    --outgroup \
    --parallel 4 --batch 10 \
    --run common
```

* On Hackintosh (i7-6700k 32G SSD)
    * `--parallel 4`: 5m57s
    * `--parallel 8`: 6m24s

* On server (E5-2690 v3 128G SSD)
    * `--parallel 4`: 6m2s
    * `--parallel 8`: 1m56s
