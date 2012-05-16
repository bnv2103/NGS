#!/bin/bash
#$ -cwd

# Call SNPs and short INDELs for one diploid individual: 
infile=$1
genome=$2
if [[ $genome == "mouse" ]];
    then
   reference="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Mus_musculus/UCSC/mm9/Sequence/BowtieIndex/genome.fa"
fi

if [[ $genome == "human" ]];
    then
    reference="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Homo_sapiens/NCBI/build37.2/Sequence/BowtieIndex/genome.fa"
fi

if [[ $genome == "rat" ]];
    then
    reference="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Rattus_norvegicus/Rattus_norvegicus/UCSC/rn4/Sequence/BowtieIndex/genome.fa"
fi

cd $infile
bamFile="./accepted_hits.bam"
echo $bamFile
echo $reference

# samtools mpileup -ugf $reference $bamFile | bcftools view -bvcg - > var.raw.bcf
# bcftools view var.raw.bcf | vcfutils.pl varFilter -D 100 > var.flt.vcf 

samtools mpileup -ugf $reference -q 5 -Q 17 $bamFile  |  bcftools view -bvcg  - > var_raw.bcf
bcftools view var_raw.bcf | vcfutils.pl varFilter -d 5 > var_flt.vcf


echo "created vcf"

#convert vcf to txt
 convert2annovar.pl var_flt.vcf -format vcf4 > var_flt_vcf
if [[ $genome == "mouse" ]];
    then
    annotate_variation.pl -geneanno -buildver mm9 var_flt_vcf /ifs/data/c2b2/ngs_lab/ngs/usr/src/annovar/mousedb
fi

if [[ $genome == "human" ]];
    then
    annotate_variation.pl -geneanno -buildver hg19 var_flt_vcf /ifs/data/c2b2/ngs_lab/ngs/usr/src/annovar/humandb
fi

if [[ $genome == "rat" ]];
    then
    annotate_variation.pl -geneanno -buildver rn4 var_flt_vcf /ifs/data/c2b2/ngs_lab/ngs/usr/src/annovar/ratdb
fi


# annotate_variation.pl -geneanno -buildver hg19 var_flt_vcf /ifs/data/c2b2/ngs_lab/ngs/usr/src/annovar/humandb

ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/mergeVPC.rb

mkdir SNPs
mv var* SNPs

#two output files are generated 
# exonic_variant_function and variant_function

# sift filtering
# annotate_variation.pl -filter -dbtype avsift -buildver hg19 var.flt.vcf.txt /ifs/data/c2b2/ngs_lab/ngs/usr/src/annovar/humandb

# SNP filtering
# annotate_variation.pl -filter -dbtype snp135 -buildver hg19 var.flt.vcf.txt /ifs/data/c2b2/ngs_lab/ngs/usr/src/annovar/humandb

# Region-based
# most conversed geomic regions

# vcfFile="var_flt_vcf"
# annotate_variation.pl -regionanno -dbtype mce46way -buildver hg19 $vcfFile /ifs/data/c2b2/ngs_lab/ngs/usr/src/annovar/humandb

# TFBs
#annotate_variation.pl -regionanno -dbtype tfbs -buildver hg19 $vcfFile /ifs/data/c2b2/ngs_lab/ngs/usr/src/annovar/humandb

# band
# annotate_variation.pl -regionanno -dbtype band -buildver hg19 $vcfFile /ifs/data/c2b2/ngs_lab/ngs/usr/src/annovar/humandb


# segmental duplications
# annotate_variation.pl -dbtype segdup -regionanno -buildver hg19 $vcfFile /ifs/data/c2b2/ngs_lab/ngs/usr/src/annovar/humandb

# previous reported structual variants in DGV
# annotate_variation.pl -regionanno -dbtype dgv -buildver hg19 $vcfFile /ifs/data/c2b2/ngs_lab/ngs/usr/src/annovar/humandb
