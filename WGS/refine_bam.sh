#!/bin/bash
#$ -cwd
	
SAMTOOLS="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/samtools"
REF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/human_g1k_v37.fasta"
chr=$1
DIR=$2"/"


$SAMTOOLS rmdup $DIR$chr.sorted.bam $DIR$chr.sorted.bam.noDup.bam
echo "Remove PCR Duplicates complete"
# rm $DIR$chr.sorted.bam
## Should add -C int : Coefficient to cap mapping quality of poorly mapped reads. See the pileup command for details. [0]
$SAMTOOLS calmd -rEAb $DIR$chr.sorted.bam.noDup.bam $REF >  $DIR$chr.sorted.bam.noDup.bam.baq.bam
$SAMTOOLS index $DIR$chr.sorted.bam.noDup.bam.baq.bam
echo "Index Complete on bam.noDup.bam.baq.bam "
rm $DIR$chr.sorted.bam.noDup.bam
echo "Realigned using Samtools- CALMD with extended BAQ done"

