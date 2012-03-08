#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -l h_vmem=1G,time=1::

DIR="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/"

##This script converts a  vcf file with one or more samples into a format required by annovar and then fires annovar's tool to get the variants.
#INPUT - $1 - a VCF file containing
#OUTPUT - $1.exome (the input require for annovar)
#	- two other files with the filtered and dropped variants.

$DIR/convert2annovar.pl $1  -format vcf4 -allallele > $1.annovar
$DIR/annotate_variation.pl --buildver hg19 $1.annovar  -filter -dbtype ljb_pp2 /ifs/data/c2b2/ngs_lab/ngs/resources/annovar_hg19/

cut -f3,4,6,7 $1.annovar.hg19_ljb_pp2_dropped > $1.annovar.a
cut -f2 $1.annovar.hg19_ljb_pp2_dropped > $1.annovar.b
echo -e "Chr\tPos\tREF\tALT\tPP2_score" > $1.pp2
paste $1.annovar.a $1.annovar.b >> $1.pp2 
rm $1.annovar.a
rm $1.annovar.b

