#!/bin/sh
#$ -cwd
## 1G,8::

prefix=$1

#samtools mpileup -b bam.full.list -uIf bcm_hg18.fasta | bcftools view -gvc - | vcfutils.pl varFilter -D 200
#samtools mpileup -b bam.full.list -uIf bcm_hg18.fasta | bcftools view -gvc -
####samtools mpileup sim1.sorted.bam -uIf bcm_hg18.fasta -C 1 -d 300 -E -q 10 -Q 15 -gDS | bcftools view -gv -e -p 0.9 -t 0.001 - 
samtools mpileup ${prefix}.sorted.bam -uIf bcm_hg18.fasta -C 1 -d 300 -E -q 10 -Q 15 -gDS > ${prefix}.bcf
bcftools view -gv -e -p 1.1 -t 0.001 ${prefix}.bcf > ${prefix}.temp.vcf
grep ^# ${prefix}.temp.vcf > ${prefix}.vcf.1.1
grep -v ^# ${prefix}.temp.vcf | awk '{if($4!="N"&&length($5)==1) print }' >> ${prefix}.vcf.1.1
rm ${prefix}.temp.vcf
cp ${prefix}.vcf.1.1 ${prefix}.vcf
#samtools mpileup simulated_reads_3.sorted.bam -uIf bcm_hg18.fasta -C 50 -d 300 -B -gDS | bcftools view -gv -e -p 1.1 -t 0.001 - 

