#!/bin/sh
#$ -cwd

len=$1
prefix=$2

# NOTE: samtools version (0.1.8 (r613)) to support "pileup" option to generate glf files
# /ifs/home/c2b2/ip_lab/yshen/usr/bin/samtools

samtools view -bS ${prefix}.sam -o ${prefix}.bam
samtools sort ${prefix}.bam ${prefix}.sorted
samtools pileup -g -T 1 -f bcm_hg18.fasta ${prefix}.sorted.bam > ${prefix}.glf
GPT_Freq -b ${prefix}.gpt.out -p 0.9 ${prefix}.glf
thunder_glf_freq --shotgun ${prefix}.gpt.out.chr21 --detailedInput -r 100 --states 200 --dosage --phase --interim 25 -o ${prefix}.final.out

