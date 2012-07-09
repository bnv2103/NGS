#!/bin/sh
#$ -cwd
## 10G,4::

chr=$1

gmap_build -d $chr -k 15 --basesize=15 ${chr}/${chr}.fasta

## --format=goby --goby-output=Sample01 Sample01.fq 

