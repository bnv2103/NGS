#!/bin/sh
#$ -cwd

i=$1

sample=`cat Sample_RK${i}.list`
goby 4g concatenate-alignments $sample -o Sample_RK${i}
goby 1g compact-file-stats Sample_RK${i}.entries

