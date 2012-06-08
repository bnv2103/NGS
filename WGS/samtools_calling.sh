#!/bin/bash
#$ -cwd

# INP="/ifs/scratch/c2b2/ac_lab/pn2204/net/bam/AC2/refine/1.noDup.rl.sorted.fxmate.bam"
INP=$1
OUT=$2
LOGS=$3

samtools mpileup -6 -ugf /ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/human_g1k_v37.fasta $INP |bcftools view -bcvg -> $OUT.bcf
# bcftools view 1.bcf | perl vcfutils.pl varFilter -D 100000000 > 1.vcf
bcftools view $OUT.bcf > $OUT.vcf

#GATK - variant annotator
echo " qsub -o $LOGS/gatk_annotate.$OUT.o -e $LOGS/gatk_annotate.$OUT.e /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/gatk_annotate_variants.sh $OUT.vcf " 
qsub -o $LOGS/gatk_annotate.$OUT.o -e $LOGS/gatk_annotate.$OUT.e /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/gatk_annotate_variants.sh $OUT.vcf 

