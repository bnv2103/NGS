#!/bin/bash
#$ -cwd
#$ -l mem=5G,time=5::

USAGE="Usage: sh _.sh -i input_bam -o output  [-m mem ]\n"
PICARD="/ifs/data/c2b2/ngs_lab/ngs/usr/picard-tools-1.65/"

while getopts m:i:p:o:g:s:h opt
  do     
  case "$opt" in
      i) input="$OPTARG";;
      p) paired="$OPTARG";;
      o) output="$OPTARG";;
      m) mem="$OPTARG";;
      g) organi="$OPTARG";;
      s) sample="$OPTARG";;
      h) echo $USAGE
	  exit 1;;
  esac
done

if [[ $organi == "" ]]
then
	organism="human"	
	echo "Assuming Human and continuing" 
else
	organism=`echo $organi | tr '[:upper:]' '[:lower:]'`  
fi
if [[ $input == "" ]]; then
    echo $USAGE
    exit 1
fi
if [[ $mem == "" ]]; then
	HEAP=2
else
	HEAP=$mem
fi
if [[ $sample == "" ]];then
	sample=`basename $input | cut -f6 -d'_' `
fi
if [[ ! $paired == "" && -e $paired ]]; then
	paired=" FASTQ2="$paired
fi
if [[ $output == "" ]]; then
	output=$input"_QC"
fi
if [[ $organism == "human" ]];then
	REF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/human_g1k_v37.fasta"
elif [[ $organism == "mouse" ]];then
	REF="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa"
elif [[ $organism == "staphylococcus" ]];then
        REF="/ifs/data/c2b2/ngs_lab/ngs/resources/genomes/Staphylococcus_aureus_USA300_FPR3757/Staphylococcus_aureus_FPR3757.fasta"
elif [[ $organism == "xenograft" ]];then
	echo "Tool not functional for $organism yet ! "
	exit
elif [[ $organism == "rat" ]]; then
	REF="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Rattus_norvegicus/Rattus_norvegicus/UCSC/rn4/Sequence/BowtieIndex/genome.fa"
elif [[ $organism == "fruitfly" ]];  then
	REF="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Drosophila_melanogaster/UCSC/dm3/Sequence/BowtieIndex/genome.fa"
elif [[ $organism == "yeast" ]];  then
	REF="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Saccharomyces_cerevisiae/Ensembl/EF3/Sequence/BowtieIndex/genome.fa"
else
        echo "Tool not functional for $organism yet ! "
        exit
fi

outputbk=$output
mkdir $output

output=$output"/$sample.picard.bam"


TEMP=$output"_temp"
if [ ! -d $TEMP ]; then
  mkdir $TEMP
fi

if [ ! -s $output ];then
	cmd="java -Xmx${HEAP}G -Djava.io.tmpdir=${TEMP} -jar ${PICARD}/FastqToSam.jar FASTQ=$input $paired  QUALITY_FORMAT=Standard OUTPUT=$output SAMPLE_NAME=$sample SORT_ORDER=coordinate "
	echo $cmd
	$cmd
fi

input=$output
output=$outputbk

cmd="java -Xmx${HEAP}G -Djava.io.tmpdir=${TEMP} -jar ${PICARD}/EstimateLibraryComplexity.jar INPUT=$input OUTPUT=$output/LibComplexity "
echo $cmd
$cmd

cmd="java -Xmx${HEAP}G -Djava.io.tmpdir=${TEMP} -jar ${PICARD}/MeanQualityByCycle.jar  INPUT=$input  OUTPUT=$output/MeanQualityByCycle REFERENCE_SEQUENCE=$REF CHART_OUTPUT=$output/MeanQualityByCycle.pdf "
echo $cmd
$cmd

cmd="java -Xmx${HEAP}G -Djava.io.tmpdir=${TEMP} -jar ${PICARD}/QualityScoreDistribution.jar INPUT=$input  OUTPUT=$output/QualityScoreDist REFERENCE_SEQUENCE=$REF CHART_OUTPUT=$output/QualityScoreDistr.pdf "
echo $cmd
$cmd

rm -rf $TEMP
