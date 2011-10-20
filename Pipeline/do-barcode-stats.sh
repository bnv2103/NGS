#!/bin/bash
#$ -cwd

files=$*

for f in $files
do 

ruby /ifs/home/c2b2/ngs_lab/ngs/code/NGS/Pipeline/barcode_stats.rb $f > $f.barcode-stats

done
