#!/bin/bash
#$ -cwd

### workflow:
#  map reads using tophat (via bowtie)
#  estimate FPKM using cufflinks with reference genes

USAGE="Usage: $0 -i fastq -s setting [-n num_threads]"

nt=4 # default 4 threads for tophat
outdir=""

while getopts i:s:n:o:h opt
  do
  case "$opt" in
      i) fq="$OPTARG";;
      s) setting="$OPTARG";;
      n) nt="$OPTARG";;
      o) outdir="$OPTARG";;
      h) echo $USAGE; exit 1
  esac
done

if [[ $fq == "" || $setting == "" ]]
    then
    echo $USAGE
    exit 1
fi

if [[ $outdir == "" ]]; then
    outdir=$fq"_tophat"
fi

. $setting

## first do tophat

# get sample name
g=`basename $fq | sed 's/.fastq//' | sed 's/\//\-/g'`
rg=$g
sampleName=$rg
library=$rg

cmd="tophat --phred64-quals --num-threads $nt  -o $outdir  --rg-id $rg --rg-sample $sampleName --rg-library $library $BOWTIEDB $fq"

echo -e "do tophat: \n $cmd"
$cmd

bam=$outdir"/accepted_hits.bam"

if [[ ! -s $bam ]]; then
    echo "tophat failed"
    exit 1
fi

## index
samtools index $bam

## do cufflinks
cuffout=$bam"_cufflinks_ref"
cuffout2=$bam"_cufflinks_ref-guided"
cmd="cufflinks -o $cuffout --GTF  $GENES --num-threads $nt  $bam"

echo -e "do cufflinks with ref genes: \n $cmd"
$cmd

cmd="cufflinks -o $cuffout2 --GTF-guide  $GENES --num-threads $nt  $bam"
echo -e "do cufflinks with ref genes as guide: \n $cmd"
$cmd

