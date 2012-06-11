#!/bin/sh
#$ -cwd

prefix=sim4

# One-time calling
# bwa index -a bwtsw bcm_hg18.fasta

# Defaults are: -b 3 -q 5 -r 2
# Expt range: -b [4-6] -q [2-4] -r [1-2]

bwa bwasw -b 4 -q 4 -r 2 -w 500 bcm_hg18.fasta ${prefix}.fq  > ${prefix}.samd
samtools calmd -S ${prefix}.samd bcm_hg18.fasta > ${prefix}.sam

samtools view -bS ${prefix}.sam chr21 | samtools sort - ${prefix}.sorted
#samtools sort ${prefix}.bam ${prefix}.sorted
samtools index ${prefix}.sorted.bam
samtools view -h ${prefix}.sorted.bam chr21 > ${prefix}.sorted.sam

#samtools pileup -g -T 1 -f bcm_hg18.fasta ${prefix}.sorted.bam > ${prefix}.glf
