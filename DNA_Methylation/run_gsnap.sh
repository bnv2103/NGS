#!/bin/sh
#$ -cwd
## 28G,1::

id=$1
chr=$2

gsnap -d $chr -k 15 -B 5 -J 0 -j 33 --mode=cmet-stranded --format=goby --goby-output=${id} ${id}.compact-reads

#gsnap -d hg19 -k 15 -B 5 --quality-protocol=sanger --mode=cmet-stranded --format=goby --goby-output=${id} ${id}.compact-reads

#gsnap -d hg19 -k 15 -B 5 --mode=cmet Sample01.fq

## --format=goby --goby-output=Sample01 Sample01.fq 

