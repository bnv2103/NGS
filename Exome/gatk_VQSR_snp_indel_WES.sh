#!/bin/bash
#$ -cwd
uname -a

heap=4
## This GATK tool requires Rscript(version 2.15.0) and ggplot2 installed and other packages dependencies. Make sure your PATH and Rlib are set as below in the ~/.bashrc and ~/.bash_profile files.
## export PATH=/ifs/data/c2b2/ngs_lab/ngs/usr/local/bin/R/R-2.15.0/bin/:/nfs/apps/R/2.14.0/bin/:$PATH
## R_LIBS=/ifs/data/c2b2/ngs_lab/ngs/usr/local/bin/R/R-2.15.0/lib64/R/library:/ifs/data/c2b2/ngs_lab/ngs/usr/local/bin/R/R-2.14.0/lib64/R/library
## R_LIBS_USER=$R_LIBS
## export R_LIBS_USER R_LIBS

while getopts v:g:m:b:h:A: opt
  do  
  case "$opt" in
      v) vcf="$OPTARG";;
      g) GLOBAL="$OPTARG";;
      m) MEM="$OPTARG";;
      b) bam="$OPTARG";;
      A) AUTO="$OPTARG";;
      h) echo $USAGE
	  exit 1;;
  esac
done

if [[ $vcf == "" || $GLOBAL == "" ]]
then
        echo $USAGE
        exit 1
fi

. $GLOBAL

if [[ ! $MEM == "" ]]; then
    heap=$MEM
fi
echo "heap size: ${heap}g"

# targets=$vcf".targets.list"
# awk '{print $1":"$2"-"$3}' $ExonFile > $targets

sh $UTILS/gatk_seperate-SNV-indel.sh -I $vcf -g $GLOBAL 
## ${RUBY18} ${UTILS}/vcf_seperate-SNV-indel.rb $vcf
 
sample_count1=`head -300 $vcf | grep "#" | tail -1 | tr "\t" "\n" | wc -l `
sample_count=$((sample_count1-9))
ADJUSTMENT=" "
INBREED=" "

if [[ $sample_count < 30 ]];then	## Set appropriate parameter for whole-exome when #samples<30
	ADJUSTMENT=" --maxGaussians 4 -percentBad 0.05 "
fi
if [[ $sample_count > 9 ]];then
	INBREED=" -an InbreedingCoeff "
fi

TEMP=$vcf"_VQSR"
mkdir -p $TEMP
JAVA="java -Xmx${heap}g -Djava.io.tmpdir="${TEMP}
GATK="$JAVA -jar "${GATKJAR}

echo $GATK
$GATK \
   -T VariantRecalibrator \
   -R $REF \
   -input $vcf.snv \
   -resource:HapMapV3,known=false,training=true,truth=true,prior=15.0 $HamMap_Sites \
   -resource:OMNI,known=false,training=true,truth=false,prior=12.0 $OMNI_Sites \
   -resource:dbSNP135,known=true,training=false,truth=false,prior=6.0 $DBSNP135 \
   -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an MQ  -an FS -an DP $INBREED $ADJUSTMENT \
   -mode SNP \
   -recalFile $vcf.SNP.recal \
   -tranchesFile $vcf.SNP.tranches \
   -rscriptFile $vcf.SNP.plots.R 
  
echo "VQSR: SNP: Variant Recalibration done."
# rm $TEMP/*

if [[ ! -s $vcf.SNP.recal  ]]; then
$GATK \
   -T ApplyRecalibration \
   -R $REF  \
   -input $vcf.snv \
   -mode SNP \
   --ts_filter_level 99.0 \
   -tranchesFile $vcf.SNP.tranches \
   -recalFile  $vcf.SNP.recal \
   -o $vcf.filered.vcf.snv

echo "VQSR: SNP: Apply Recalibration done."
fi

rm -rf $TEMP

## separate indel from SNV:
# ${RUBY18} ${UTILS}/vcf_seperate-SNV-indel.rb $vcf.filtered.vcf
sh ${BPATH}/do_release.sh $vcf.filtered.vcf.snv

## filter indels:
dname=`dirname $vcf`
if [[ $AUTO == "" ]]; then
	cmd="qsub -l mem=6G,time=24:: -N filter-indel -o $dname/filter-indel.o -e $dname/filter-indel.e ${BPATH}/vcf_filter-indel.sh -I $vcf.indel -s $vcf.filered.vcf.snv -g $GLOBAL "
else
        cmd="qsub -l mem=6G,time=24:: -N filter-indel -o $dname/filter-indel.o -e $dname/filter-indel.e ${BPATH}/vcf_filter-indel.sh -I $vcf.indel -s $vcf.filered.vcf.snv -g $GLOBAL  -A AUTO"
fi
echo $cmd
$cmd

