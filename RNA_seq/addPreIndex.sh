#!/bin/bash
#$ -cwd

infile=$1

pre=$2

# cp $isoforms $outdir"/"$infileB".bam"

mv $1 $pre"_"$1