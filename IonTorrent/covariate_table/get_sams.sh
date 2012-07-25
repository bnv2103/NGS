#!/bin/bash
# -cwd

# directory containing bam files
bamdir=`pwd`/$1
# target list
bed=$2

for s in `ls -d $bamdir/*.bam`
do
    outfile=`pwd`/`basename $s | cut -f1 -d "."`.sam
    echo samtools view -q 30 -L $bed $s | qsub -l mem=2G,time=1:: -o $outfile -e `pwd`/get_sams.e
done
