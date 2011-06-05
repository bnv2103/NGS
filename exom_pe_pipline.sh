#!/bin/bash
### Author : H Ryan Kim
### Version : 0.7
### Year : 2010

###### Flow Chart ###############################
### 1. Raw Fastq
### 2. Convert illumina Fastq to Sanger Fastq
### 3. QC (FastQC)
### 4. Alignment (BWA)
### 5. Convert SAM to BAM (Picard)
### 6. Local Realignment around Indels (GATK)
### 7. Remove Duplicates
### 8. Base Quality Score Recalibration
### 9. Call SNPs and InDels (GATK)
#################################################

###### Prerequite
### 1. Raw Fastq (convert qseq to Fastq or convert Colorspace to Letterspace(for SOLiD))
### 2. Convert illumina Fastq to Sanger Fastq (covert illumina quality value to sanger's)
### 3. QC (FastQC) (fastqc - generate qc report)


###### Usage:
##  sh exom_pe_pipline.sh 1stFastq 2ndFastq ProjectName Reference ROD

###### User Defined Inputs
Fastq_1=$1 
Fastq_2=$2
Project=$3
REF=$4
ROD=$5

###### Environmental Variables
TMPDIR=$PWD/tmp
GATK=/ifs/data/shares/deepsequencing/solid/src/GenomeAnalysisTK-1.0.4168/GenomeAnalysisTK.jar
GATKacov=/ifs/data/shares/deepsequencing/solid/src/GenomeAnalysisTK-1.0.4168/AnalyzeCovariates.jar
GATKR=/ifs/data/shares/deepsequencing/solid/src/GATK/Sting/R
Rbin=/nfs/apps/R/2.9.0/bin/Rscript
GATKPERL=/ifs/data/shares/deepsequencing/solid/src/GATK/Sting/perl
GATKPYTHON=/ifs/data/shares/deepsequencing/solid/src/GATK/Sting/python
PICARD=/ifs/data/shares/deepsequencing/solid/src/picard-tools-1.43/
mkdir $PWD/tmp

#### 4. Alignment (BWA)
## Align with BWA
bwa aln -t 8 $REF $Fastq_1 > $Fastq_1.sai;
bwa aln -t 8 $REF $Fastq_2 > $Fastq_2.sai;
## Generate alignment in the SAM format
bwa sampe $REF $Fastq_1.sai $Fastq_2.sai $Fastq_1 $Fastq_2 > $Project.sam; 


#### 5. Convert SAM to BAM (Picard)
## Sort bwa SAM file using PICARD TOOLS SortSam.jar - this will also produce the BAM file
java -jar  SortSam.jar SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT \ 
INPUT=$Project.sam OUTPUT=$Project.bam;
## samtools index
samtools index $Project.bam


#### 6. Local Realignment around Indels (GATK)
##  6.1. Creating Intervals : RealignerTargetCreator 
java –Xmx5g -jar $GATK -T RealignerTargetCreator -R $REF -D $ROD \
-I  $Project.bam -o $Project.bam.forRealigner.intervals;

##  6.2. Realigning : IndelRealigner 
java -Djava.io.tmpdir=$TMPDIR –Xmx5g -jar $GATK -T IndelRealigner \
-R $REF -D $ROD -I  $Project.bam -o $Project.realn.bam \
-targetIntervals $Project.bam.forRealigner.intervals;
## samtools index
samtools index $Project.realn.bam;

## 6.3. Sort realigned BAM file using PICARD TOOLS SortSam.jar 
## GATK IndelRealigner  produces a name sorted BAM
java –Xmx6g -jar $PICARD/SortSam.jar \
INPUT=$Project.realn.bam OUTPUT=$Project.realn.sorted.bam \
SORT_ORDER=coordinate TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=SILENT;
## samtools index
samtools index $Project.realn.soretd.bam;

## 6.4. Fixing the mate pairs of realigned reads using Picard tools FixMateInformation.jar
java -Djava.io.tmpdir=$TMPDIR -jar -Xmx6g FixMateInformation.jar \
INPUT=$Project.realn.sorted.bam OUTPUT=$Project.realn.sorted.fixed.bam \
SO=coordinate VALIDATION_STRINGENCY=SILENT TMP_DIR=$TMPDIR;
## samtools index
samtools index $Project.realn.sorted.fixed.bam ;


### 7. Remove Duplicates
## Remove duplicate reads with Picard tools MarkDuplicates.jar 
java -Xmx6g –jar MarkDuplicates.jar \
INPUT=$Project.realn.sorted.fixed.bam \
OUTPUT=$Project.realn.duperemoved.bam \
METRICS_FILE=$Project.realn.Duplicate.metrics.file \
REMOVE_DUPLICATES=true \
ASSUME_SORTED=false TMP_DIR=$TMPDIR \
VALIDATION_STRINGENCY=SILENT;
## samtools index
samtools index $Project.realn.duperemoved.bam;


### 8. Base Quality Score Recalibration
## 8.1. GATK CountCovariates 
java -Xmx8g -jar $GATK -T CountCovariates -R $REF --DBSNP $ROD \
-I $Project.realn.duperemoved.bam \
-recalFile $Project.realn.duperemoved.bam.recal_data.csv \
--default_platform Illumina \
-cov ReadGroupCovariate \
-cov QualityScoreCovariate \
-cov CycleCovariate \
-cov DinucCovariate \
-cov TileCovariate \
-cov HomopolymerCovariate \
-nback 5;

## 8.2. GATK AnalyzeCovariates
java -Xmx5g –jar $GATKacov \
-recalFile $Project.realn.duperemoved.bam.recal_data.csv \
-outputDir $Project.realn.duperemoved.bam.recal.plots \
-resources $GATKR  \
-Rscript $Rbin;

## 8.3. GATK TableRecalibration 
java –Xmx6g -jar $GATK -T TableRecalibration -R $REF \
-I $Project.realn.duperemoved.bam \
--out $Project.final.bam \
-recalFile $Project.realn.duperemoved.bam.recal_data.csv \
--default_platform Illumina;

##samtools index
samtools index $Project.final.bam;


### 9. Call SNPs and InDels (GATK)

## SNP Calling 
java -Xmx5g -jar $GATK -T UnifiedGenotyper -R $REF -D $ROD \
-baq CALCULATE_AS_NECESSARY -baqGOP 30 -nt 8 \
-A DepthOfCoverage -A AlleleBalance -A HaplotypeScore -A HomopolymerRun -A MappingQualityZero -A QualByDepth -A RMSMappingQuality -A SpanningDeletions \
-I $Project.final.bam -o $Project.raw.snps.vcf \
-verbose $Project.raw.snps.vcf.verbose -metrics $Project.raw.snps.vcf.metrics;
## VariantFiltration  & annotation
java –Xmx5g -jar $GATK -T VariantFiltration -R $REF -D $ROD \
-o $Project.VariantFiltration.snps.vcf \
-B variant,VCF, $Project.raw.snps.vcf \
-B mask,Bed, indels.mask.bed --maskName InDel \
--clusterSize 3 --clusterWindowSize 10 \
--filterExpression "DP <= 8" --filterName "DP8" \
--filterExpression "SB > -0.10" --filterName "StrandBias" \
--filterExpression "HRun > 8" --filterName "HRun8" \
--filterExpression "QD < 5.0" --filterName "QD5" \
--filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" --filterName "hard2validate"; 

##  Generate raw indel calls
java -Xmx6g -jar $GATK -T IndelGenotyperV2 -R $REF --DBSNP $ROD \
-I $Project.final.bam \
-bed $Project.raw.indels.bed \
-o $Project.raw.indels.detailed.output.vcf  \
--metrics_file $Project.gatk.raw.indels.metrics.file \
-verbose $Project.gatk.raw.indels.verbose.output.bed \
-minCoverage 8 -S SILENT –mnr 1000000;
## Filter raw indels
perl  $GATKPERL/filterSingleSampleCalls.pl --calls $Project.gatk.raw.indels.verbose.output.bed \
--max_cons_av_mm 3.0 --max_cons_nqs_av_mm 0.5 --mode ANNOTATE > $Project.gatk.filtered.indels.bed


