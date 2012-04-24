#!/bin/bash
#$ -cwd

# print basic stats
out=$1
fileName=$out".csv"

ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/printTitle_PE.rb $fileName
for f in *cufflinks; do ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/comb_stats_PE.rb $f $fileName; done
