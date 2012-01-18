#!/bin/bash
#$ -cwd

# print basic stats
out=$1

 fileName=$out".pdf"
 Rscript /ifs/scratch/c2b2/ngs_lab/xs2182/code/printMA.R $fileName
 echo "done MA"
