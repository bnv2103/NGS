#!/bin/bash
#$ -cwd

genome=$1

for f in *1.fastq*; do
    f3=`echo $f | sed 's/_1.fastq/_3.fastq/' `
    if [ -e $f3 ];then sh /ifs/scratch/c2b2/ngs_lab/xs2182/code/fusionGene.sh $genome $f $f3; fi

    # if [ ! -e $f3 ]; then
    # 	sh /ifs/scratch/c2b2/ngs_lab/xs2182/code/RNA_pipeline.sh $genome $f;
    # fi

done