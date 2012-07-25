#!/bin/bash
#$ -cwd

USAGE="Usage: sh _.sh -i input_bam -o output  [-m mem ]\n"
PICARD="/ifs/data/c2b2/ngs_lab/ngs/usr/picard-tools-1.65/"

WIN=200	#Window size = 2x(template size)  (default 200 for HiSeq ) for Miseq =300

## Collects GC bias , InsertSize metrics from an Aligned/MApped BAM file.
## Reference file is required in the setting's file as REF="blah.fasta"
## Window size is optional.
while getopts m:i:o:g:w:h opt
  do     
  case "$opt" in
      i) input="$OPTARG";;
      o) output="$OPTARG";;
      m) mem="$OPTARG";;
      g) setting="$OPTARG";;
      w) WIN="$OPTARG";;
      h) echo $USAGE
	  exit 1;;
  esac
done

if [[ ! $setting == "" ]];then
	. $setting
fi
if [[ $input == "" ]]; then
    echo $USAGE
    exit 1
fi

if [[ $output == "" ]]; then
	output=$input
fi

if [[ $mem == "" ]]; then
	HEAP=2
else
	HEAP=$mem
fi

TEMP=$input"_GCtemp"
if [ ! -d $TEMP ]; then
  mkdir $TEMP
fi
if [ ! -d $output ];then
	mkdir $output
fi

cmd="java -Xmx${HEAP}G -Djava.io.tmpdir=${TEMP} -jar ${PICARD}/CollectGcBiasMetrics.jar INPUT=$input  OUTPUT=$output.GCbias_detail CHART=$output.GCbias.pdf REFERENCE_SEQUENCE=$REF VALIDATION_STRINGENCY=SILENT WINDOW_SIZE=$WIN "
echo $cmd
$cmd

cmd="java -Xmx${HEAP}G -Djava.io.tmpdir=${TEMP} -jar ${PICARD}/CollectInsertSizeMetrics.jar INPUT=$input OUTPUT=$output.InsertSize_detail HISTOGRAM_FILE=$output.InsertSize.pdf REFERENCE_SEQUENCE=$REF  VALIDATION_STRINGENCY=SILENT  "
echo $cmd
$cmd

rm -rf $TEMP
