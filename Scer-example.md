# Build alignments of several strains of *Saccharomyces cerevisiae*

As example for egaz and alignDB. Extract from Scer_wgs of [`OPs-download.md`](https://github.com/wang-q/withncbi/blob/master/pop/OPs-download.md) and [`OPs-align.md`](https://github.com/wang-q/withncbi/blob/master/pop/OPs-align.md).

## Download

1. Create `scer_example.tsv` manually.

    ```bash
    mkdir -p ~/data/alignment/example/scer      # operation directory
    mkdir -p ~/data/alignment/example/GENOMES   # sequence directory

    cd ~/data/alignment/example/GENOMES

    echo -e '#name\tprefix\torganism\tcontigs' > example_wgs.tsv
    echo -e "Spar\tAABY\tSaccharomyces paradoxus NRRL Y-17217\t832" >> example_wgs.tsv
    ```

2. Create working directory and download WGS sequences.

    ```bash
    cd ~/data/alignment/example/GENOMES

    perl ~/Scripts/withncbi/taxon/wgs_prep.pl \
        -f example_wgs.tsv \
        --fix \
        -o WGS \
        -a

    aria2c -UWget -x 6 -s 3 -c -i WGS/example_wgs.url.txt

    find WGS -name "*.gz" | xargs gzip -t
    ```

3. Download strains of *Saccharomyces cerevisiae* at good assembly status.

    ```bash
    mkdir -p ~/data/alignment/example/GENOMES/DOWNLOAD
    cd ~/data/alignment/example/GENOMES/DOWNLOAD

    # Download S288c, EC1118 and RM11_1a separately
    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000146045.2.assembly.txt \
        --nuclear -name S288c \
        > S288c.seq.csv

    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCA_000218975.1.assembly.txt \
        --nuclear --genbank --scaffold -name EC1118 \
        > EC1118.seq.csv

    perl ~/Scripts/withncbi/taxon/assembly_csv.pl \
        -f ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCA_000149365.1.assembly.txt \
        --genbank --scaffold -name RM11_1a \
        > RM11_1a.seq.csv

    echo "#strain_name,accession,strain_taxon_id,seq_name" > example.seq.csv
    cat S288c.seq.csv EC1118.seq.csv RM11_1a.seq.csv \
        | perl -nl -e '/^#/ and next; /^\s*$/ and next; print;' \
        >> example.seq.csv

    # Download, rename files and change fasta headers
    perl ~/Scripts/withncbi/taxon/batch_get_seq.pl \
        -p -f example.seq.csv
    ```

## Align

In these steps, I don't use ensembl annotations. Just Assembly and RepeatMasker ones.

Stats of pairwise and multiple alignments are minimal.

1. `gen_pop_conf.pl`

    ```bash
    mkdir -p ~/data/alignment/example/scer
    cd ~/data/alignment/example/scer

    perl ~/Scripts/withncbi/pop/gen_pop_conf.pl \
        -i ~/data/alignment/example/GENOMES/WGS/example_wgs.data.yml \
        -o scer_test.yml \
        -d ~/data/alignment/example/GENOMES/WGS \
        -m prefix \
        -r '*.fsa_nt.gz' \
        --opt group_name=scer \
        --opt base_dir='~/data/alignment/example' \
        --opt data_dir="~/data/alignment/example/scer" \
        --opt rm_species=Fungi \
        --dd ~/data/alignment/example/GENOMES/DOWNLOAD \
        --download "name=S288c;taxon=559292" \
        --download "name=RM11_1a;taxon=285006" \
        --download "name=EC1118;taxon=643680" \
        --plan 'name=Scer_n3_Spar;t=S288c;qs=EC1118,RM11_1a,Spar;o=Spar' \
        --plan 'name=Scer_n3_pop;t=S288c;qs=EC1118,RM11_1a' \
        --plan 'name=Scer_n2_Spar;t=S288c;qs=RM11_1a,Spar;o=Spar' \
        -y
    ```

2. Rest routing things.

    ```bash
    cd ~/data/alignment/example/scer

    # pop_prep.pl
    perl ~/Scripts/withncbi/pop/pop_prep.pl -p 12 -i scer_test.yml

    sh 01_file.sh
    sh 02_rm.sh
    sh 03_strain_info.sh

    # plan_ALL.sh
    sh plan_ALL.sh

    sh 1_real_chr.sh
    sh 3_pair_cmd.sh
    sh 4_rawphylo.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    # other plans
    sh plan_Scer_n3_Spar.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    sh plan_Scer_n3_pop.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh

    sh plan_Scer_n2_Spar.sh
    sh 5_multi_cmd.sh
    sh 7_multi_db_only.sh
    ```

## Copy directories and files to example directory

```bash
cd ~/Scripts/alignDB

cp -R ~/data/alignment/example/scer/Genomes/S288c data/

cp -R ~/data/alignment/example/scer/Pairwise/S288cvsRM11_1a data/
cp -R ~/data/alignment/example/scer/Pairwise/S288cvsEC1118 data/
cp -R ~/data/alignment/example/scer/Pairwise/S288cvsSpar data/

cp ~/data/alignment/example/scer/id2name.csv data/
cp ~/data/alignment/example/scer/fake_tree.nwk data/
```

## Build local ensembl database

Use `build_ensembl.pl`.

If don't use scripts in `gene/*`, local ensembl database names are arbitrary.

```bash
cd ~/Scripts/alignDB

cp -R ~/data/ensembl82/mysql/saccharomyces_cerevisiae_core_29_82_4 data/

perl util/build_ensembl.pl --checksum --ensembl data/saccharomyces_cerevisiae_core_29_82_4
perl util/build_ensembl.pl --initdb --db yeast --ensembl data/saccharomyces_cerevisiae_core_29_82_4
perl util/build_ensembl.pl --initdb --db saccharomyces_cerevisiae_core_29_82_4 --ensembl data/saccharomyces_cerevisiae_core_29_82_4
```

## Two-way alignments

`-taxon` and `-chr` are omitted because they already exist in default ones.

`--ensembl yeast` means different in step 20 and 22. In step 20, `yeast` is the mysql database name.
And in step 22, `yeast` is an alias to `saccharomyces_cerevisiae_core_29_82_4`.

```bash
cd ~/Scripts/alignDB/data

# S288cvsRM11_1a
perl ~/Scripts/alignDB/extra/two_way_batch.pl \
    -t "559292,S288c" -q "285006,RM11_1a" \
    -d S288cvsRM11_1a \
    -da S288cvsRM11_1a \
    -lt 5000 \
    --parallel 8 \
    -r all

# S288cvsSpar
perl ~/Scripts/alignDB/extra/two_way_batch.pl \
    -t "559292,S288c" -q "226125,Spar" \
    -d S288cvsSpar \
    --ensembl yeast \
    -da S288cvsSpar \
    -lt 5000 \
    --parallel 8 \
    -r all
```

## Three-way alignments by `join_dbs.pl`

```bash
cd ~/Scripts/alignDB/data

# without outgroup
perl ~/Scripts/alignDB/extra/join_dbs.pl \
    --no_insert --block --trimmed_fasta --length 1000 \
    --goal_db S288cvsRM11_1avsSpar \
    --target 0target \
    --queries 0query,1query \
    --dbs S288cvsRM11_1a,S288cvsSpar

# with outgroup
perl ~/Scripts/alignDB/extra/join_dbs.pl \
    --no_insert --block --trimmed_fasta --length 1000 \
    --goal_db S288cvsRM11_1arefSpar \
    --target 0target \
    --queries 0query \
    --outgroup 1query \
    --dbs S288cvsRM11_1a,S288cvsSpar

find . -type f -name "*.fas" | parallel -j 8 gzip
```

## Three-way alignments by multiz

```bash
cd ~/Scripts/alignDB/data

if [ -d ScervsRM11_1a_Spar_mz ]; then
    rm -fr ScervsRM11_1a_Spar_mz;
    mkdir -p ScervsRM11_1a_Spar_mz;
else
    mkdir -p ScervsRM11_1a_Spar_mz;
fi;

if [ -d ScervsRM11_1a_Spar_fasta ]; then
    rm -fr ScervsRM11_1a_Spar_fasta;
fi;

if [ -d ScervsRM11_1a_Spar_refined ]; then
    rm -fr ScervsRM11_1a_Spar_refined;
fi;

# mz
perl ~/Scripts/egaz/mz.pl \
    -d S288cvsRM11_1a \
    -d S288cvsSpar \
    --tree fake_tree.nwk \
    --out ScervsRM11_1a_Spar_mz \
    -p 8
find ScervsRM11_1a_Spar_mz -type f -name "*.maf" | parallel --no-run-if-empty -j 8 gzip

# maf2fas
mkdir -p ScervsRM11_1a_Spar_fasta
find ScervsRM11_1a_Spar_mz -name "*.maf" -or -name "*.maf.gz" \
    | parallel --no-run-if-empty -j 8 fasops maf2fas {} -o ScervsRM11_1a_Spar_fasta/{/}.fas

# refine fasta
perl ~/Scripts/egaz/refine_fasta.pl \
    --msa mafft --block -p 8 \
    --outgroup \
    -i ScervsRM11_1a_Spar_fasta \
    -o ScervsRM11_1a_Spar_refined

find ScervsRM11_1a_Spar_refined -type f -name "*.fas" | parallel -j 8 gzip

rm -fr ScervsRM11_1a_Spar_mz
rm -fr ScervsRM11_1a_Spar_fasta
```

## Multi-way batch

```bash
cd ~/Scripts/alignDB/data

perl ~/Scripts/alignDB/extra/multi_way_batch.pl \
    -d ScervsRM11_1a_Spar \
    -da ScervsRM11_1a_Spar_refined \
    --ensembl yeast \
    --block \
    --id id2name.csv \
    --outgroup \
    -lt 1000 --parallel 8 --batch 5 \
    --run all
```

## Slicing

### Yeast intergenic regions

* two-way (axt)

```bash
mkdir -p ~/Scripts/alignDB/data/feature
cd ~/Scripts/alignDB/data/feature
 
perl ~/Scripts/alignDB/slice/write_runlist_feature.pl \
    -d S288cvsRM11_1a -e yeast --feature intergenic -l 500

mkdir axt
mv S288cvsRM11_1a.intergenic.yml axt
rm axt/*.axt

perl ~/Scripts/alignDB/slice/write_align_slice.pl \
    -d S288cvsRM11_1a -f axt/S288cvsRM11_1a.intergenic.yml -t axt

perl ~/Scripts/alignDB/extra/two_way_batch.pl \
    -d S288cvsRM11_1a_intergenic \
    -t "559292,S288c" -q "285006,RM11_1a" \
    -da axt \
    -lt 1000 \
    --parallel 8 \
    --run basic
```

* multi-way (fas)

```bash
mkdir -p ~/Scripts/alignDB/data/feature
cd ~/Scripts/alignDB/data/feature
 
perl ~/Scripts/alignDB/slice/write_runlist_feature.pl \
    -d ScervsRM11_1a_Spar -e yeast --feature intergenic -l 500

mkdir fas
mv ScervsRM11_1a_Spar.intergenic.yml fas
rm fas/*.fas

perl ~/Scripts/alignDB/slice/write_align_slice.pl \
    -d ScervsRM11_1a_Spar -f fas/ScervsRM11_1a_Spar.intergenic.yml -t fas --outgroup

perl ~/Scripts/alignDB/extra/multi_way_batch.pl \
    -d S288cvsRM11_1a_intergenic \
    --id ../id2name.csv \
    -da fas \
    --block \
    -lt 1000 \
    --parallel 8 \
    --run basic
```

### Yeast non-repeat regions

```bash
mkdir -p ~/Scripts/alignDB/data/feature
cd ~/Scripts/alignDB/data/feature
 
perl ~/Scripts/alignDB/slice/write_runlist_feature.pl \
    -d S288cvsRM11_1a -e yeast --feature repeat -l 500

```
