#!/bin/bash
#$ -cwd

bam=$1
outfile=$2
# sizeFile=$3
# chr=$4
# chrSeq=$5

chr="1"
sizeFile="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Homo_sapiens/NCBI/build37.2/Sequence/BowtieIndex/chr.sizes"
chrSeq="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Homo_sapiens/NCBI/build37.2/Sequence/BowtieIndex/1.fa"

 samtools sort accepted_hits.bam accepted_hits.sorted
 samtools index accepted_hits.sorted.bam

java -Xmx2g -jar /ifs/data/c2b2/ngs_lab/ngs/usr/src/scripture-beta2.jar -alignment accepted_hits.sorted.bam -out $outfile -sizeFile $sizeFile -chr $chr -chrSequence $chrSeq  

# scripture-beta2.jar