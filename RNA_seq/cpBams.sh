#!/bin/bash
#$ -cwd

infile=$1
outdir=$2

mkdir $outdir
# copy bams generated from cufflink fold to current working directory
 isoforms=$infile"/accepted_hits.bam"
# isoforms=$infile"/accepted_hits.sorted.bam"
# isoforms1=$infile"/accepted_hits.sorted.bam.bai"

infileB=` echo $infile | cut -f6 -d '_'`
cp $isoforms $outdir"/"$infileB".bam"
# cp $isoforms1 $outdir"/"$infileB".bam.bai"