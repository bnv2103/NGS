#!/bin/sh
#$ -cwd
### 18G,1::

sample=$1
chr=$2

goby 30g discover-sequence-variants  variants/Sample_RK${sample} --format genotypes --groups A=Sample_RK${sample} -o genotypes${sample}.vcf --genome $chr

echo "METHYLATION COMPLETE"

