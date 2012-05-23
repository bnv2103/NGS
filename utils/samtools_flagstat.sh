#!/bin/bash
#$ -cwd

bam=$1
samtools flagstat $bam > $bam.flagstat
