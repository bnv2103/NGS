#!/bin/bash
#$ -cwd

geno=$1 #moues or human ...
group1=$2 #sample1,sample2,sample3
group2=$3 #sample4,sample5,sample6
dir=$4 # working directory
lab1=$5 # label 1 sample1_sample2_sample3
lab2=$6 # label 2 sample4_samle5_sample6

cd $dir
output=$group1"_"$group2
mkdir $output
ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/pipeline_cuffdiff.rb $geno $group1 $group2 $dir $output

mv $output*.txt $output

if [[ $geno == "mouse" ]];
    then
    genome="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Mus_musculus/UCSC/mm9/Sequence/BowtieIndex/genome.fa"
    gtf="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf"
fi

if [[ $geno == "human" ]];
    then
    genome="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex/genome.fa"     
    gtf="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
fi

cd $output
cuffmerge -g $gtf -s $genome -p 1 *transcripts.txt
bams=`cat *bams.txt`
echo $bams

cuffdiff -o "cuffdiff/" -L $lab1,$lab2 -u merged_asm/merged.gtf $bams


cd cuffdiff
Rscript /ifs/scratch/c2b2/ngs_lab/xs2182/code/getSig.R gene_exp.diff $output
Rscript /ifs/scratch/c2b2/ngs_lab/xs2182/code/getSig.R isoform_exp.diff $output

cp $output* ../../


