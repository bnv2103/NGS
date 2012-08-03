#!/bin/sh
#$ -cwd
### mem=6G,time=8::

# One time calling only

base=$1
fileA=$2
fileB=$3

goby 5g alignment-to-annotation-counts ${base}*.entries --annotation mm10-exons-ensembl-annotation.txt --include-annotation-types gene --compare A/B --groups A=`cat ${fileA} | tr '\n' ','`/B=`cat ${fileB} | tr '\n' ','` --stats stats.tsv
