#!/bin/bash
#$ -cwd

### workflow:
#  map reads using tophat (via bowtie)
#  estimate FPKM using cufflinks with reference genes

USAGE="Usage: $0 -f fastqF -r fastqR -s setting [-n num_threads]"

nt=4 # default 4 threads for tophat
outdir=""

while getopts f:r:s:n:o:h opt
  do
  case "$opt" in
      f) fq1="$OPTARG";;
      r) fq2="$OPTARG";;
      s) setting="$OPTARG";;
      n) nt="$OPTARG";;
      o) outdir="$OPTARG";;
      h) echo $USAGE; exit 1
  esac
done

if [[ $fq1 == "" || $fq2 == "" || $setting == "" ]]
    then
    echo $USAGE
    exit 1
fi

if [[ $outdir == "" ]]; then
    outdir=$fq1"_tophat"
fi

. $setting

## first do tophat

# get sample name
g=`basename $fq1 | sed 's/.fastq//' | sed 's/\//\-/g'`
rg=$g
sampleName=$rg
library=$rg

cmd="tophat --phred64-quals --num-threads $nt  -o $outdir  --rg-id $rg --rg-sample $sampleName --rg-library $library --mate-inner-dist 200 $BOWTIEDB $fq1 $fq2"

echo -e "do tophat: \n $cmd"
$cmd

bam=$outdir"/accepted_hits.bam"

if [[ ! -s $bam ]]; then
    echo "tophat failed"
    exit 1
fi


## do cufflinks
cuffout=$bam"_cufflinks_ref-guide"

cmd="cufflinks -o $cuffout --compatible-hits-norm --GTF-guide  $GENES $bam"

echo -e "do cufflinks without ref genes: \n $cmd"
$cmd

## do cufflinks with guide
cuffout=$bam"_cufflinks_ref"

cmd="cufflinks -o $cuffout --compatible-hits-norm --GTF  $GENES $bam"

echo -e "do cufflinks with ref genes: \n $cmd"
$cmd

# after run cufflink, statistics of number of reads and FPKM
ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/comb_stats.rb $outdir "summary.csv"
