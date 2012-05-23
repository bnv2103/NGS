#!/bin/bash
#$ -cwd

infile=$1
outfile=$2


ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/countLines_fastq.rb $infile $outfile