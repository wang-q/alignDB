#!/bin/bash
cd ~/data/yeast10000

# prepare reference
# cat ~/data/alignment/yeast65/S288C/*.fa > ref/S288C.fa.masked
# perl -p -e '/>/ and next; $_ = uc' ref/S288C.fa.masked > ref/S288C.fa
# rm ref/S288C.fa.masked

# index reference genome
# bwa index -a bwtsw ref/S288C.fa
# samtools faidx ref/S288C.fa

# prepare vcf
# perl ~/Scripts/alignDB/gvf2vcf.pl ~/data/ensembl65/gvf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.gvf > S288C.unsort.vcf
# ~/share/vcftools/vcf-sort S288C.unsort.vcf > ref/S288C.vcf

# rice reference
# cat ~/data/alignment/rice/nip_65/chr[1-9].fa ~/data/alignment/rice/nip_65/chr1[0-2].fa> ref/nip_65.fa.masked
# perl -p -e '/>/ and next; $_ = uc' ref/nip_65.fa.masked  > ref/nip_65.fa
# rm ref/nip_65.fa.masked
# ~/bin/x86_64/faSize  -detailed ref/nip_65.fa > ref/chr.sizes

# QC
# find . -name "*.sam.gz" | xargs fastqc -f sam -t 6

# sra to fastq (pair end)
~/share/sratoolkit/fastq-dump ERP000547/ERR038793/ERR038793.sra \
    --split-files --gzip -O process/ERR038793 

# align each samples to genome
bwa aln -t 4 -q 15 ref/S288C.fa process/ERR038793/ERR038793_1.fastq.gz \
    > process/ERR038793/ERR038793_1.sai
bwa aln -t 4 -q 15 ref/S288C.fa process/ERR038793/ERR038793_2.fastq.gz \
    > process/ERR038793/ERR038793_2.sai

# convert sai to sam
# add read groups info
bwa sampe -r "@RG\tID:ERR038793\tLB:ERR038793\tPL:ILLUMINA\tSM:ERR038793" \
    ref/S288C.fa \
    process/ERR038793/ERR038793*.sai process/ERR038793/ERR038793*.fastq.gz \
    > process/ERR038793/aln.sam

# conver sam to bam
samtools view -bS process/ERR038793/aln.sam -o process/ERR038793/aln.bam

# sort bam
samtools sort process/ERR038793/aln.bam process/ERR038793/aln.sorted

# index bam
samtools index process/ERR038793/aln.sorted.bam

# merge with samtools
# samtools merge finalBamFile.bam *.bam

# merge with picard
# java  ¨CXmx1g -jar ~/share/picard/MergeSamFiles.jar \
#    INPUT=SM_IDa.sorted  INPUT=SM_IDb.sorted  INPUT=SM_IDn.sorted  \
#    USE_THREADING=true  VALIDATION_STRINGENCY=LENIENT \
#    AS=true OUTPUT=SM.sorted.bam 

#java -Xmx2g -jar /gpfssan1/home/wangq/share/imrdenom/external//MergeSamFiles.jar \
#    TMP_DIR=/gpfssan1/home/wangq/data/yeast10000/process_imr/ERR038793/A/  \
#    I=/gpfssan1/home/wangq/data/yeast10000/process_imr/ERR038793/A/aln_A1.sp1.bam \
#    O=/gpfssan1/home/wangq/data/yeast10000/process_imr/ERR038793/A/mapset_withdup_0.bam \
#    ASSUME_SORTED=true \
#    VALIDATION_STRINGENCY=SILENT \
#    SO=coordinate \
#    > /gpfssan1/home/wangq/data/yeast10000/process_imr/ERR038793/A/rs_mg.0 2>&1

# index regions for realignment
java -Xmx1g -jar ~/share/GenomeAnalysisTK/GenomeAnalysisTK.jar -nt 4 \
    -T RealignerTargetCreator \
    -R ref/S288C.fa \
    -I process/ERR038793/aln.sorted.bam  \
    --out process/ERR038793/aln.intervals

# realign bam to get better Indel calling
java -Xmx1g -jar ~/share/GenomeAnalysisTK/GenomeAnalysisTK.jar \
    -T IndelRealigner \
    -R ref/S288C.fa \
    -I process/ERR038793/aln.sorted.bam \
    -targetIntervals process/ERR038793/aln.intervals \
    --out process/ERR038793/aln.realigned.bam

# dup marking
java -Xmx1g -jar ~/share/picard/MarkDuplicates.jar \
    INPUT=process/ERR038793/aln.realigned.bam \
    OUTPUT=process/ERR038793/aln.dedup.bam \
    METRICS_FILE=process/ERR038793/output.metrics \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=LENIENT

# reindex the realigned dedup BAM
samtools index process/ERR038793/aln.dedup.bam

# recalibration - Count covariates
java -Xmx1g -jar ~/share/GenomeAnalysisTK/GenomeAnalysisTK.jar -nt 4 \
    -T CountCovariates  \
    -R ref/S288C.fa \
    -I process/ERR038793/aln.dedup.bam \
    -knownSites ref/S288C.vcf \
    -recalFile process/ERR038793/recal_data.csv \
    -cov ReadGroupCovariate \
    -cov QualityScoreCovariate \
    -cov CycleCovariate \
    -cov DinucCovariate 
    
# recalibration - Tabulate recalibration
java -Xmx1g -jar ~/share/GenomeAnalysisTK/GenomeAnalysisTK.jar \
    -T TableRecalibration  \
    -R ref/S288C.fa \
    -I process/ERR038793/aln.dedup.bam \
    -o process/ERR038793/aln.recal.bam \
    -recalFile process/ERR038793/recal_data.csv

# recalibration - Analyze Covariates
# java -Xmx1g -jar ~/share/GenomeAnalysisTK/AnalyzeCovariates.jar \
#   -recalFile process/ERR038793/recal_data.csv  \
#   -outputDir process/ERR038793/  \
#   -ignoreQ 5

# reindex the realigned BAM
# gatk automaticly generates a bai
# samtools index process/ERR038793/aln.recal.bam

# coverage
#java -Xmx1g -jar ~/share/GenomeAnalysisTK/GenomeAnalysisTK.jar \
#    -T DepthOfCoverage \
#    -R ref/S288C.fa \
#    -I process/ERR038793/aln.recal.bam \
#    --minMappingQuality 20 \
#    -o process/ERR038793/coverage.txt 

## dindel stage 1
#[% bin_dir.dindel %]/dindel --analysis getCIGARindels \
#    --bamFile [% item.dir %]/[% item.name %].baq.bam \
#    --ref [% ref_file.seq %] \
#    --outputFile [% item.dir %]/[% item.name %].dindel_output
#
## dindel stage 2
#[% bin_dir.dindel %]/makeWindows.py \
#    --inputVarFile [% item.dir %]/[% item.name %].dindel_output.variants.txt \
#    --windowFilePrefix [% item.dir %]/[% item.name %].realign_windows \
#    --numWindowsPerFile 1000
#
## dindel stage 3
#[% bin_dir.dindel %]/dindel --analysis indels --doDiploid \
#    --bamFile [% item.dir %]/[% item.name %].baq.bam \
#    --ref [% ref_file.seq %] \
#    --varFile [% item.dir %]/[% item.name %].realign_windows.2.txt \
#    --libFile [% item.dir %]/[% item.name %].dindel_output.libraries.txt \
#    --outputFile [% item.dir %]/[% item.name %].dindel_stage2_output_windows.2
#
## dindel stage 4
#[% bin_dir.dindel %]/mergeOutputDiploid.py \
#    --inputFiles [% item.dir %]/[% item.name %].dindel_stage2_outputfiles.txt \
#    --ref [% ref_file.seq %] \
#    --outputFile [% item.dir %]/[% item.name %].variantCalls.VCF

# generate fastq from bam
samtools mpileup -uf ref/S288C.fa process/ERR038793/aln.recal.bam \
    | bcftools view -cg - \
    | vcfutils.pl vcf2fq > process/ERR038793/aln.fq

# convert fastq to fasta 
# mask bases with quality lower than 20 to lowercases
seqtk fq2fa process/ERR038793/aln.fq 20 > process/ERR038793/aln.fa

# let's clean up
find process/ERR038793/ -type f \
    -name "*.sai"    -o -name "*.fastq.gz" \
    -o -name "*.bam" -o -name "*.bai" \
    -o -name "*.sam" -o -name "*.intervals" \
    -o -name "*.csv" -o -name "*.metrics" \
    | grep -v "recal" | xargs rm

bsub -n 8 -J velopt -q largemem ~/share/VelvetOptimiser/VelvetOptimiser.pl -s 27 -e 50 -t 4 -f ' -shortPaired -bam MiSeq_Ecoli_MG1655_110721_PF.bam '