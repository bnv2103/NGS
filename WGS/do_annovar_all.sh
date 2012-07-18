#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -l h_vmem=8G,time=2::

DIR="/ifs/scratch/c2b2/ngs_lab/sz2317/softwares/annovar/"

##This script converts a  vcf file with one or more samples into a format required by annovar and then fires annovar's tool to get the variants.
#INPUT - $1 - a VCF file containing
#OUTPUT - $1.annovar (the input require for annovar)
#	- two other files with the filtered and dropped variants.
echo "Begin Annovar Preprocessing : Input file $1";
if [ ! -e $1.annovar ] ;then
	$DIR/convert2annovar.pl $1  -format vcf4 -includeinfo -allallele > $1.annovar
	echo "Converted VCF to Annovar compatible file : $1.annovar";
fi

echo "Begin Annovar Annotation";
$DIR/summarize_annovar1.pl  $1.annovar $DIR"/humandb/" -outfile $1.summary -ver1000g 1000g2012feb -verdbsnp 132 -genetype=refgene --buildver hg19

echo "Finished Annovar Annotation : $1.summary.*";

# $DIR/annotate_variation.pl --buildver hg19 $1.annovar  -filter -dbtype ljb_pp2 /ifs/data/c2b2/ngs_lab/ngs/resources/annovar_hg19/
# $DIR/annotate_variation.pl --buildver hg19 -filter -dbtype avsift $1.annovar  /ifs/data/c2b2/ngs_lab/ngs/resources/annovar_hg19/ &

# Convert annovar result back to VCF using original VCF for header
echo "Converting Annovar Summary to VCF : $1.summary.exome/genome_summary.csv";
perl /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/GIT_CODE/NGS/WGS/convert_annovar_vcf-all-samples.pl $1.summary.exome_summary.csv $1
perl /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/GIT_CODE/NGS/WGS/convert_annovar_vcf-all-samples.pl $1.summary.genome_summary.csv $1
echo "Annovar COMPLETE ";
echo "OUTPUT - $1.summary.exome_summary.csv.vcf ";
echo "OUTPUT - $1.summary.genome_summary.csv.vcf ";

