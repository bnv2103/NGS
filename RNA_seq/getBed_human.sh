#!/bin/bash
#$ -cwd

infile=$1
outdir=$2
# cd $infile
# bamFile="accepted_hits.bam"

 for j in {1..22}
 do
    bamToBed -i $infile | awk -v j="$j" '{if ( $1 == j){print "chr"$0;} }' > $outdir"/"chr$j.bed
    
 done

  bamToBed -i $infile | awk '{if ( $1 == "Y"){print "chr"$0;} }' > $outdir"/"chrY.bed
  bamToBed -i $infile | awk '{if ( $1 == "X"){print "chr"$0;} }' > $outdir"/"chrX.bed


cd $outdir
mkdir bedFiles
mv chr*.bed bedFiles
