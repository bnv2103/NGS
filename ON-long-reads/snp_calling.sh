#!/bin/sh
#$ -cwd
## 1G,8::

prefix=$1
pval=$2
suffix=$3
region=$4

if [[ $region -eq 0 ]]
then
	region=$SGE_TASK_ID
fi
rbase=${prefix}_$region
vbase=${suffix}_$region

samtools mpileup ${rbase}.sorted.bam -uIf reference/bcm_hg18.fasta -C 1 -d 300 -E -q 10 -Q 15 -gDS > ${vbase}.bcf
bcftools view -gv -e -p $pval -t 0.001 ${vbase}.bcf > ${vbase}.temp.vcf
grep ^# ${vbase}.temp.vcf > ${vbase}.${pval}.vcf
grep -v ^# ${vbase}.temp.vcf | awk -F'\t' '{if($4!="N"&&length($5)==1) { split($8,info,";");split(info[1],dp,"="); if(dp[2]>=2) print $0;} }' >> ${vbase}.${pval}.vcf
rm ${vbase}.temp.vcf

#samtools mpileup -b bam.full.list -uIf bcm_hg18.fasta | bcftools view -gvc - | vcfutils.pl varFilter -D 200
#samtools mpileup -b bam.full.list -uIf bcm_hg18.fasta | bcftools view -gvc -
####samtools mpileup sim1.sorted.bam -uIf bcm_hg18.fasta -C 1 -d 300 -E -q 10 -Q 15 -gDS | bcftools view -gv -e -p 0.9 -t 0.001 - 

