#!/bin/bash
#$ -cwd

infile=$1
outdir=$2
# cd $infile
# bamFile="accepted_hits.bam"

 for j in {1..19}
 do
    bamToBed -i $infile | awk -v j="$j" '{if ( $1 == "chr"j){print $0;} }' > $outdir"/"chr$j.bed
    
 done
  bamToBed -i $infile | awk '{if ( $1 == "chrM"){print $0;} }' > $outdir"/"chrM.bed
  bamToBed -i $infile | awk '{if ( $1 == "chrY"){print $0;} }' > $outdir"/"chrY.bed
  bamToBed -i $infile | awk '{if ( $1 == "chrX"){print $0;} }' > $outdir"/"chrX.bed

mkdir bedFiles
mv chr*.bed bedFiles
# mkdir bedFiles
# mv *.bed bedFiles
