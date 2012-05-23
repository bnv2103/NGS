#!/bin/bash
#$ -cwd


for f in *cufflinks; do ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/generateTranscripts.rb $f; done

if [[ $1 == "mouse" ]];
    then
    genome="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Mus_musculus/UCSC/mm9/Sequence/BowtieIndex/genome.fa"
    gtf="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf"
fi

if [[ $1 == "human" ]];
    then
    genome="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Homo_sapiens/NCBI/build37.2/Sequence/BowtieIndex/genome.fa"
    gft="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Homo_sapiens/NCBI/build37.2/Annotation/Genes/genes.gtf"
fi



cuffmerge -g $gtf -s $genome -p 1 transcripts.txt
