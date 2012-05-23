#!/bin/bash
#$ -cwd

infile=$1
outdir=$2

# copy bams generated from cufflink fold to current working directory
isoforms=$infile"/accepted_hits.bam"


infileB=` echo $infile | cut -f6 -d '_'`
cp $isoforms $outdir"/"$infileB".bam"
