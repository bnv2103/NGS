#!/bin/bash
#$ -cwd

bam=$1
nt=$2

if [[ $nt == "" ]]; then
     nt=4
fi

outdir=$bam"_no-ref"
cufflinks -o $outdir --num-threads $nt  $bam

