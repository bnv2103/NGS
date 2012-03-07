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

fname="$1.annovar"

cut -f3,4,6,7 $fname.hg19_ljb_pp2_dropped > $fname.a
cut -f2 $fname.hg19_ljb_pp2_dropped > $fname.b
echo -e "Chr\tPos\tREF\tALT\tPP2_score" > $fname.pp2
paste $fname.a $fname.b >> $fname.pp2 
rm $fname.a
rm $fname.b

