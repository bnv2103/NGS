#!/bin/sh
#$ -cwd
## 8G,8::

id1=$1
id2=$2
chr=$3

goby 2g discover-sequence-variants Sample*.entries --output meth.vcf --compare A/B --groups A=Sample01,Sample02/B=Sample03,Sample04  --genome v37 --format methylation

#kgoby 3g discover-sequence-variants --output meth.vcf ${id1} ${id2} --compare A/B --groups A=${id1}/B=${id2}  --genome $chr --diploid true

##--format methylation
## --compare A/B

