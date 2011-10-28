#!/bin/bash
#$ -cwd

# print basic stats
out=$1
isMA=$2
fileName=$out".csv"

ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/printTitle.rb $fileName
for f in *fastq_tophat ; do ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/comb_stats.rb $f $fileName; done

echo "done calculation"

dirName="pdfPlot"
mkdir $dirName
for f in *fastq_tophat ; do  sh /ifs/scratch/c2b2/ngs_lab/xs2182/code/cpIsoforms.sh  $f $dirName; done
cd $dirName

fileName=$out".pdf"
Rscript /ifs/scratch/c2b2/ngs_lab/xs2182/code/printPDF.R $fileName  $isMA 
echo "done PDF"