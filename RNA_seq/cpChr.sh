#!/bin/bash
#$ -cwd

infile=$1
outdir=$2

mkdir $outdir
# copy bams generated from cufflink fold to current working directory
files=$infile"/chr6"


infileB=` echo $infile | cut -f6 -d '_'`
cp $files $outdir"/"$infileB"_chr6.txt"

files=$infile"/chr7"


infileB=` echo $infile | cut -f6 -d '_'`
cp $files $outdir"/"$infileB"_chr7.txt"

