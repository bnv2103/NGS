#!/bin/bash
#$ -cwd

# two parameters: out = output file name, isMA = 0(no MA plot), 1 (need MA plot)

out=$1
isMA=$2
echo "printing pdf"

 dirName="to_release"
 if [ ! -d $dirName ]; then mkdir -p $dirName; fi

 for f in *cufflinks; do  sh /ifs/scratch/c2b2/ngs_lab/xs2182/code/cpIsoforms.sh  $f $dirName; done
 # for f in *cufflinks; do  sh /ifs/scratch/c2b2/ngs_lab/xs2182/code/cpIsoforms_nonRef.sh  $f $dirName; done
 cd $dirName

 fileName=$out".pdf"
 Rscript /ifs/scratch/c2b2/ngs_lab/xs2182/code/printPDF.R $fileName  $isMA
 echo "done PDF"

 rm *genes
 rm *isoforms
# rm *nonRef
cp ../../summary.csv ./
cp ../summary.csv ./
cp /ifs/scratch/c2b2/ngs_lab/xs2182/code/README* ./