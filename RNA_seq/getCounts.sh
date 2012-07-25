#!/bin/bash
#$ -cwd

samFile=$1
txtFile=$2
f=$3
gtf=$4

samtools view -h -o $samFile $f

htseq-count $samFile $gtf --stranded=no > $txtFile
