#!/bin/bash
#$ -cwd

dir=$1 #input dir - where is the bam



samtools view $dir"/accepted_hits.bam" | cut -f3 | sort | uniq -c > $dir"/reads"

         
Rscript /ifs/scratch/c2b2/ngs_lab/xs2182/code/printSpikeIn.R $dir"/reads" $dir"_trans"
