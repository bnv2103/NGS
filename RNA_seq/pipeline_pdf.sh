#!/bin/bash
#$ -cwd

# two parameters: out = output file name, isMA = 0(no MA plot), 1 (need MA plot)

out=$1
isMA=$2

RNABASE="/ifs/scratch/c2b2/ngs_lab/xs2182/code"
echo "printing pdf"

 dirName="to_release"
 if [ ! -d $dirName ]; then mkdir -p $dirName; fi
 # ruby $RNABASE/printTitle.rb "summary.csv"
 # for f in *cufflinks; do ruby $RNABASE/comb_stats.rb $f "summary_test.csv"; done

 for f in *cufflinks; do  sh $RNABASE/cpIsoforms.sh  $f $dirName; done
 # for f in *cufflinks; do  sh /ifs/scratch/c2b2/ngs_lab/xs2182/code/cpIsoforms_nonRef.sh  $f $dirName; done
 cd $dirName

 Rscript $RNABASE/combineFPKM.R
 fileName=$out".pdf"
 Rscript $RNABASE/printPDF.R $fileName  $isMA
 echo "done PDF"

 rm *genes
 rm *isoforms
# rm *nonRef
cp ../../summary.csv ./
cp ../summary.csv ./
cp $RNABASE/README* ./

# for f in *cufflinks; do
# qsub -l mem=1G,time=6:: /ifs/scratch/c2b2/ngs_lab/xs2182/code/cpBams.sh $f $dirName
# done