#!/bin/bash
#$ -cwd

bam=$1
out=$2
samtools view $bam 7:141950000-142550000 > $out

