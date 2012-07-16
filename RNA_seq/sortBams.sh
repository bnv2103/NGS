#!/bin/bash
#$ -cwd

infile=$1

cd $infile
samtools sort accepted_hits.bam accepted_hits.sorted
samtools index "accepted_hits.sorted.bam"

