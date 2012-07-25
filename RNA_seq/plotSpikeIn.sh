#!/bin/bash
#$ -cwd

# print basic stats
for f in *csv;
do Rscript /ifs/scratch/c2b2/ngs_lab/xs2182/code/calcSpikeTranRatio.R $f $f
done

