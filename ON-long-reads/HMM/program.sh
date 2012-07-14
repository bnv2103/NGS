#!/bin/sh
#$ -cwd
### 30G,4::

prefix=$1
suffix=$2
region=$3
pval=$4
out=$5

if [[ $region -eq 0 ]]
then
	region=$SGE_TASK_ID
fi

rbase=${prefix}_${region}.sorted
vbase=${suffix}_${region}.${pval}.vcf
obase=${out}_${region}.${pval}

echo $rbase $vbase $region $obase

./program $rbase $vbase $region $obase

