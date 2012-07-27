#!/bin/bash
#$ -cwd

samFile=$1
txtFile=$2
bamFile=$3
gtf=$4

samtools view -h -o $samFile $bamFile

htseq-count $samFile $gtf --stranded=no > $txtFile
