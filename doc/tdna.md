# T-DNA insertion sites

## Download

Wed Jan 20 15:03:28 2016

```bash
mkdir -p ~/data/salk
cd ~/data/salk

# Alternative url: http://signal.salk.edu/database/
perl ~/Scripts/download/list.pl -u http://natural.salk.edu/database/ -r transcriptome
perl ~/Scripts/download/download.pl -i database.yml -a
```

## Convert to beds

```bash
mkdir -p ~/data/salk/Atha
mkdir -p ~/data/salk/OstaJap

for name in CMT CSHL FLAG GABI IMAL MX RATM SAIL SALK SK WISC
do
    echo $name;
    perl ~/Scripts/alignDB/ofg/tdna2bed.pl -i ~/data/salk/database/tdnaexpress/T-DNA.$name -o ~/data/salk/Atha/T-DNA.$name.bed;
done

for name in affjp cirad csiro genoplante gsnu ostid PFG_FSTs rifgp rmd ship trim ucd
do
    echo $name;
    perl ~/Scripts/alignDB/ofg/tdna2bed.pl -i ~/data/salk/database/RiceGE/T-DNA.$name -o ~/data/salk/OstaJap/T-DNA.$name.bed;
done

```

## Restore databases

```bash
perl ~/Scripts/alignDB/util/dup_db.pl -g Athavsself_tdna -f ~/data/dumps/mysql/Athavsself.sql.gz

```

## Use style 'center'

```bash
perl ~/Scripts/alignDB/ofg/insert_bed.pl \
    -d Athavsself_tdna  --style center --batch 10 --parallel 8 \
    --tag tdna --type CMT  -f ~/data/salk/Atha/T-DNA.CMT.bed  \
    --tag tdna --type CSHL -f ~/data/salk/Atha/T-DNA.CSHL.bed \
    --tag tdna --type FLAG -f ~/data/salk/Atha/T-DNA.FLAG.bed \
    --tag tdna --type GABI -f ~/data/salk/Atha/T-DNA.GABI.bed \
    --tag tdna --type IMAL -f ~/data/salk/Atha/T-DNA.IMAL.bed \
    --tag tdna --type MX   -f ~/data/salk/Atha/T-DNA.MX.bed   \
    --tag tdna --type RATM -f ~/data/salk/Atha/T-DNA.RATM.bed \
    --tag tdna --type SAIL -f ~/data/salk/Atha/T-DNA.SAIL.bed \
    --tag tdna --type SALK -f ~/data/salk/Atha/T-DNA.SALK.bed \
    --tag tdna --type SK   -f ~/data/salk/Atha/T-DNA.SK.bed   \
    --tag tdna --type WISC -f ~/data/salk/Atha/T-DNA.WISC.bed

perl ~/Scripts/alignDB/init/update_sw_cv.pl -d Athavsself_tdna --batch 10 --parallel 8
perl ~/Scripts/alignDB/init/update_feature.pl -d Athavsself_tdna -e ath_65 --batch 10 --parallel 8

perl ~/Scripts/alignDB/stat/ofg_stat_factory.pl --by tt -d Athavsself_tdna -o Athavsself_tdna.ofg.xlsx

```
