#!/bin/sh
#$ -cwd
## 6G,:15:

chr=$1

goby 5g build-sequence-cache --basename $chr ${chr}/${chr}.fasta.gz

## --format=goby --goby-output=Sample01 Sample01.fq 

