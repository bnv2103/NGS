#!/bin/sh
#$ -cwd

#samtools mpileup -b bam.full.list -uIf bcm_hg18.fasta | bcftools view -gvc - | vcfutils.pl varFilter -D 200
#samtools mpileup -b bam.full.list -uIf bcm_hg18.fasta | bcftools view -gvc -
samtools mpileup simulated_reads_4.sorted.bam -uIf bcm_hg18.fasta -C 1 -d 300 -E -q 10 -Q 15 -gDS | bcftools view -gv -e -p 1.1 -t 0.001 - 
#samtools mpileup simulated_reads_3.sorted.bam -uIf bcm_hg18.fasta -C 50 -d 300 -B -gDS | bcftools view -gv -e -p 1.1 -t 0.001 - 

