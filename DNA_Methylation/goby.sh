#!/bin/sh
#$ -cwd
## 8G,8::

id1=$1
id2=$2
chr=$3

# Note that supplying file names that include a folder path can be problematic here given the separation of groups using a "/".
# Alternative is to copy over all aligned files (.entries, .header etc) from all folders to a common one and then run.

#goby 4g discover-sequence-variants --output meth.vcf Sample01 Sample02 --compare A/B --groups A=Sample01/B=Sample02  --genome v37 --diploid true --format methylation

goby 3g discover-sequence-variants --output meth.vcf ${id1} ${id2} --compare A/B --groups A=${id1}/B=${id2}  --genome $chr --diploid true

##--format methylation
## --compare A/B

