#!/bin/sh
#$ -cwd
### 18G,1::

file=$1
file1=$2
file2=$3
chr=$4

merged_file="merged.list"

goby 30g discover-sequence-variants `cat ${merged_file}` --compare A/B --groups A=`cat ${file1} | tr '\n' ','`/B=`cat ${file2} | tr '\n' ','`  --format methylation --output meth.vcf --genome $chr

echo "METHYLATION COMPLETE"

