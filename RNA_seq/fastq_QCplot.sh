#!/bin/bash
#$ -cwd
 
infile=$1 # input fastq file: sample.fastq
outfile=$2 # output pdf file : sample.pdf

ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/fastq_QCplot.rb $infile $outfile


