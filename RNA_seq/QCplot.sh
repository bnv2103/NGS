#!/bin/bash
#$ -cwd

infileDir=$1 # input fastq file dir : infileDir for sample.bam
# outfile=$2 # output pdf file : infileDir/sample.pdf

# cd $infileDir
ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/QCplot.rb $infileDir 


