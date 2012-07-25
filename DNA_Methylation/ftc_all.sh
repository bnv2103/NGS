#!/bin/sh
#$ -cwd
## 2G,8::

dir=$1

cd $1
goby 8g ftc --quality-encoding Sanger `cat file.list`
cd ..

#goby 1g ftc --quality-encoding Sanger --verbose-quality-scores Sample01.fq Sample02.fq 

## --genome hg19 --compare A/B --format methylation --output meth.vcf

