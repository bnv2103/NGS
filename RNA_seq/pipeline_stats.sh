#!/bin/bash
#$ -cwd

# print basic stats, create csv files for isoforms and genes, print histogram and MA plots. 
# isMA = 1, if need MA plots 
out=$1
isMA=$2
fileName=$out".csv"

ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/printTitle.rb $fileName
for f in *cufflinks; do ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/comb_stats.rb $f $fileName; done

echo "done calculation"

dirName="result_to_release"
mkdir $dirName
for f in *cufflinks; do  sh /ifs/scratch/c2b2/ngs_lab/xs2182/code/cpIsoforms.sh  $f $dirName; done
cd $dirName
echo "copied individual cufflink results" 

fileName=$out".pdf"
Rscript /ifs/scratch/c2b2/ngs_lab/xs2182/code/printPDF.R $fileName  $isMA 
echo "done PDF"

cp ../$out".csv" ./

rm *genes
rm *isoforms
cp ../../summary.csv ./
cp ../summary.csv ./
cp /ifs/scratch/c2b2/ngs_lab/xs2182/code/README* ./
