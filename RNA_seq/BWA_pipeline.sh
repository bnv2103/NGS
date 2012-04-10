#!/bin/bash
#$ -cwd

# REF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/human.fasta"
 REF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/human_g1k_v37.fasta"
fastq1=$1

maxgaps=2
maxeditdist=0.04
qualtrim=5
platform="illumina"
threads=4

# bwa aln -q $qualtrim -o $maxgaps -n $maxeditdist -t  $threads  $REF  $fastq1 > $fastq1".sai"

bwa aln -t $threads $REF $fastq1 > $fastq1".sai"

echo "done sai"

bwa samse $REF $fastq1".sai" $fastq1 > $fastq1".sam"

echo "done sam"

# cmd1="bwa samse $REF $fastq1.sai $fastq1 | $samtools view -bS  >  $fastq1.bam"




  