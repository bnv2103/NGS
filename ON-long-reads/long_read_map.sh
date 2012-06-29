#!/bin/sh
#$ -cwd
### 6G,8::

prefix=$1

# One-time calling
# bwa index -a bwtsw bcm_hg18.fasta

# Defaults are: -b 3 -q 5 -r 2
# Expt range: -b [4-6] -q [2-4] -r [1-2]

bwa bwasw -b 4 -q 4 -r 2 -w 500 bcm_hg18.fasta ${prefix}.fq  > ${prefix}.samd
samtools calmd -S ${prefix}.samd bcm_hg18.fasta > ${prefix}.sam

grep -v "NN" ${prefix}.sam | awk 'function rem(cigar) {num=0;for(it=1;it<=length(cigar);it++) {sb=substr(cigar,it,1); if(sb ~ /[0-9]/) {num=10*num+sb;} else return num;}} {if((NF==3)||($4>0&&$3=="chr21"&&rem($6)<2000)) print}' > temp.sam
mv ${prefix}.sam ${prefix}_temp.sam
mv temp.sam ${prefix}.sam

samtools view -bS ${prefix}.sam | samtools sort - ${prefix}.sorted
samtools index ${prefix}.sorted.bam
samtools view -h ${prefix}.sorted.bam > ${prefix}.sorted.sam

