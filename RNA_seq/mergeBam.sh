#!/bin/bash
#$ -cwd

infile1=$1
infile2=$2

samtools merge merged.bam $infile1 $infile2