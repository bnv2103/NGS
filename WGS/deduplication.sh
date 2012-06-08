#!/bin/bash
#$ -cwd

HEAP=4

while getopts i:s:m:d:h o
  do 
  case "$o" in
      i) bam="$OPTARG";;
      s) setting="$OPTARG";;
      m) mem="$OPTARG";;
      d) dir="$OPTARG";;
      h) echo $USAGE
	  exit 1;;
  esac
done


if [[ $bam == "" ||  $setting == "" ]]; then
    echo "Usage: $0 -i foo.bam -s setting [ -m heap ]"
    exit 0
fi

.  $setting

if [[ ! $MEM == "" ]]; then
    HEAP=$MEM
fi

temp=$dir"/"$bam".temp_dir/"
basebam=`basename $bam`
outbam="$dir/$basebam"
mkdir -p $temp

JAVA="java -Xmx${HEAP}g -Djava.io.tmpdir="${temp}
MarkDup="$JAVA -jar ${PICARD}/MarkDuplicates.jar"


$MarkDup I=$bam O="$outbam.dedup.bam" METRICS_FILE="$outbam.dup" VALIDATION_STRINGENCY=SILENT

if [[ $? == 0 ]]
    then
#    mv $outbam.dedup.bam $outbam
#    rm $bam.bai
    $SAMTOOLS index "$outbam.dedup.bam"
    
    echo "dedup complete"
    rm -rf $temp
else
    echo "dedup failed"
    rm -rf $temp
    exit 1
fi
