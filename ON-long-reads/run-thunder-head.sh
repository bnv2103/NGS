#!/bin/sh
#$ -cwd

len=$1
prefix=$2

head -${len} ${prefix}.sam > ${prefix}_head.sam
samtools view -bS ${prefix}_head.sam -o ${prefix}_head.bam
samtools sort ${prefix}_head.bam ${prefix}_head.sorted
samtools pileup -g -T 1 -f bcm_hg18.fasta ${prefix}_head.sorted.bam > ${prefix}_head.glf
GPT_Freq -b ${prefix}_head.gpt.out -p 0.9 ${prefix}_head.glf
thunder_glf_freq --shotgun ${prefix}_head.gpt.out.chr21 --detailedInput -r 100 --states 200 --dosage --phase --interim 25 -o ${prefix}_head.final.out

