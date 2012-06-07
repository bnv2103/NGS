#!/bin/bash
#$ -cwd

infile=$1

# infileB=` echo $infile | cut -f3,4 -d '_'`

# cp $isoforms $outdir"/"$infileB".bam"

cp $1 $2/$1