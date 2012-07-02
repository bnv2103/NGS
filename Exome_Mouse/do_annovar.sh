#!/bin/bash
#$ -cwd

ANNOVAR="/ifs/scratch/c2b2/ngs_lab/sz2317/softwares/annovar/"
BPATH="/ifs/scratch/c2b2/ngs_lab/ngs/code/NGS/Exome_Mouse/"

if [ ! -e $1.annovar ]
then
	perl $ANNOVAR/convert2annovar.pl $1  -format vcf4 -allallele --includeinfo > $1.annovar
fi
 
perl $ANNOVAR/annotate_variation.pl -filter -dbtype snp128 --buildver mm9  $1.annovar  $ANNOVAR/mousedb/ -webfrom annovar 

cat $1.mm9_snp128_dropped | awk 'BEGIN{OFS="\t";}{print $3,$4,$5,$6,$7,$8,$9,$2,$11,$12,$13,$14,$15";"$1,$16,$17; }' > $1.dbsnp.annovar
cat $1.mm9_snp128_filtered >> $1.dbsnp.annovar

perl $ANNOVAR/annotate_variation.pl --geneanno  -dbtype gene --buildver mm9  $1.dbsnp.annovar  $ANNOVAR/mousedb/ -webfrom annovar 

wait

perl $BPATH/annovar_geneanno_vcf.pl $1.dbsnp.annovar.variant_function $1.dbsnp.annovar.exonic_variant_function $1 > $1.complete.annotated.vcf


