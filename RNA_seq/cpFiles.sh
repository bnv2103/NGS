#!/bin/bash
#$ -cwd

infile=$1

infileB=` echo $infile | cut -f6,8 -d '_'`

# cp $isoforms $outdir"/"$infileB".bam"

mv $1 $2/$infileB