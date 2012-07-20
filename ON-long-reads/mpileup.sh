#!/bin/sh
#$ -cwd

rbase=$1
region=$2
vbase=$3
pval=$4
i=$SGE_TASK_ID

lent=`samtools view -H ${rbase}.sorted.bam | grep -w "SN:chr${region}" | cut -d':' -f 3`
let "overlap=$lent/200"

if [[ $i -ne 1 ]]
then
	let "start=(${i}-1)*$lent/10+1-$overlap"
else
	let "start=(${i}-1)*$lent/10+1"
fi

if [[ $i -ne 10 ]]
then
	let "end=${i}*$lent/10+$overlap"
else
	let "end=${i}*$lent/10"
fi

echo $i $lent $start $end

samtools mpileup ${rbase}.sorted.bam -uIf reference/bcm_hg18.fasta -C 1 -d 300 -E -q 10 -Q 15 -gDS -r "chr${region}:${start}-${end}" > ${vbase}.bcf.$i

bcftools view -gvN -e -p $pval -t 0.001 ${vbase}.bcf.${i} > ${vbase}.${i}.temp.vcf

