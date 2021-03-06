#!/bin/bash
#$ -cwd

#  gtf_mouse="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf"
#  gtf_human="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Homo_sapiens/NCBI/build37.2/Annotation/Genes/genes.gtf"
genome=$1

if [[ $genome == "mouse" ]];
    then
   gtf="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf"
fi

if [[ $genome == "human" ]];
    then
    gtf="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Homo_sapiens/NCBI/build37.2/Annotation/Genes/genes.gtf"
fi


for f in *sam

do

 samFile=$f
 txtFile=$f".txt"

 qsub -l mem=2G,time=5:: /ifs/scratch/c2b2/ngs_lab/xs2182/code/getCounts.sh $samFile $txtFile $f $gtf


done
