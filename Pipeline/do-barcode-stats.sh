#!/bin/bash
#$ -cwd

files=$*

for f in $files
do 

ruby /ifs/home/c2b2/ngs_lab/ngs/code/NGS/Pipeline/barcode_stats.rb $f > $f.barcode-stats

done

for i in `seq 1 8` ; do echo "lane $i" ; head fastq/s_"$i"n_2.fastq.barcode-stats ; done > top_barcode
