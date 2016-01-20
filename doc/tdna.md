# T-DNA insertion sites

## Download

Wed Jan 20 15:03:28 2016

```bash
mkdir -p ~/data/salk
cd ~/data/salk

# Alternative url: http://signal.salk.edu/database/
perl ~/Scripts/download/list.pl -u http://natural.salk.edu/database/
perl ~/Scripts/download/download.pl -i database.yml -a
```

## Convert to beds

```bash
mkdir -p ~/data/salk/ath
mkdir -p ~/data/salk/nip

for name in CMT CSHL FLAG GABI IMAL MX RATM SAIL SALK SK WISC
do
    echo $name;
    perl ofg/tdna2bed.pl -i ~/data/salk/database/tdnaexpress/T-DNA.$name -o ~/data/salk/ath/T-DNA.$name.bed;
done

for name in affjp cirad csiro genoplante gsnu ostid PFG_FSTs rifgp rmd ship trim ucd
do
    echo $name;
    perl ofg/tdna2bed.pl -i ~/data/salk/database/RiceGE/T-DNA.$name -o ~/data/salk/nip/T-DNA.$name.bed;
done

```

## Create vsself

