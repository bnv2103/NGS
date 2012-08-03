#!/bin/sh
#$ -cwd

i=$1

sample=`cat Sample_RK${i}.list`
goby 7g concatenate-alignments $sample -o Sample_RK${i}
goby 3g compact-file-stats Sample_RK${i}.entries

