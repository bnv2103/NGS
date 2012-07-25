#!/bin/bash
#$ -cwd

infile=$1 # input fastq file: sample.fastq


ruby /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/underWork/Pipeline/fastq_extractSeq.rb CTGTAGGCACCATCAATAGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG $infile

