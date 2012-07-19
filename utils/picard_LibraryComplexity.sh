#!/bin/bash
#$ -cwd

USAGE="Usage: sh _.sh -i input_bam -o output  [-m mem ]\n"
PICARD="/ifs/data/c2b2/ngs_lab/ngs/usr/picard-tools-1.65/"
HEAP=2

while getopts m:i:p:o:s:h opt
  do     
  case "$opt" in
      i) input="$OPTARG";;
      p) paired="$OPTARG";;
      o) output="$OPTARG";;
      m) mem="$OPTARG";;
      s) sample="$OPTARG";;
      h) echo $USAGE
	  exit 1;;
  esac
done

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
	sample=`basename $input `
fi
sample=" SAMPLE_NAME=$sample "

if [[ ! $paired == "" && -e $paired ]]; then
	paired=" FASTQ2="$paired
fi
if [[ $output == "" ]]; then
	output=$input".LibComplexity"
fi

TEMP=$output"_temp_lib_comp"
if [ ! -d $TEMP ]; then
  mkdir $TEMP
fi

if [ ! -s $output.bam ];then
	cmd="java -Xmx${HEAP}G -Djava.io.tmpdir=${TEMP} -jar ${PICARD}/FastqToSam.jar FASTQ=$input $paired  QUALITY_FORMAT=Standard OUTPUT=$output.bam $sample SORT_ORDER=coordinate "
	echo $cmd
	$cmd
fi

cmd="java -Xmx${HEAP}G -Djava.io.tmpdir=${TEMP} -jar ${PICARD}/EstimateLibraryComplexity.jar INPUT=$output.bam  OUTPUT=$output "
echo $cmd
$cmd
rm -rf $TEMP
rm $output.bam


######### Future , create own tool
# cmd="java -Xmx${HEAP}G -Djava.io.tmpdir=${TEMP} -jar ${PICARD}/MeanQualityByCycle.jar  INPUT=$input  OUTPUT=$output/MeanQualityByCycle REFERENCE_SEQUENCE=$REF CHART_OUTPUT=$output/MeanQualityByCycle.pdf "
# echo $cmd
# $cmd

# cmd="java -Xmx${HEAP}G -Djava.io.tmpdir=${TEMP} -jar ${PICARD}/QualityScoreDistribution.jar INPUT=$input  OUTPUT=$output/QualityScoreDist REFERENCE_SEQUENCE=$REF CHART_OUTPUT=$output/QualityScoreDistr.pdf "
# echo $cmd
# $cmd

