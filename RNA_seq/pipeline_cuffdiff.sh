#!/bin/bash
#$ -cwd

genome=$1
group1=$2 #sample1,sample2,sample3
group2=$3 #sample4,sample5,sample6
dir=$4

output=$group1"_"$group2
mkdir $output
ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/pipeline_cuffdiff.rb $genome $group1 $group2 $dir $output

mv $output*.txt $output

if [[ $genome == "mouse" ]];
    then
    genome="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Mus_musculus/UCSC/mm9/Sequence/BowtieIndex/genome.fa"
    gtf="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf"
fi

if [[ $genome == "human" ]];
    then
    genome="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Homo_sapiens/UCSC/hg19/Sequence/BowtieIndex/genome.fa"     
    gtf="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
fi

cd $output
cuffmerge -g $gtf -s $genome -p 1 *transcripts.txt
bams=`cat *bams.txt`
echo $bams
cuffdiff -o "cuffdiff/" -L $group1,$group2 $gtf $bams
cd cuffdiff
Rscript /ifs/scratch/c2b2/ngs_lab/xs2182/code/getSig.R gene_exp.diff $output
Rscript /ifs/scratch/c2b2/ngs_lab/xs2182/code/getSig.R isoform_exp.diff $output

cp $output* ../../

