#!/bin/sh
#$ -cwd
## 10G,8::

id1=$1
id2=$2

goby 8g ftc --quality-encoding Sanger ${id1} ${id2}

#goby 1g ftc --quality-encoding Sanger --verbose-quality-scores Sample01.fq Sample02.fq 

## --genome hg19 --compare A/B --format methylation --output meth.vcf

