#!/bin/bash
#$ -cwd

bam=$1
setting=$2
nt=$3

if [[ $nt == "" ]]; then
     nt=4
fi


. $setting

outdir=$bam"_RABT"
cufflinks -o $outdir --num-threads $nt --GTF-guide  $GENES  $bam
