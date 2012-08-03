#!/bin/sh
#$ -cwd
### mem=2G,time=:5:

#Doesn't really require qsub

i=$1

sample=Sample_RK${i}
goby 1g coverage $sample -a mm10-capture-regions -o ${sample}.tsv

