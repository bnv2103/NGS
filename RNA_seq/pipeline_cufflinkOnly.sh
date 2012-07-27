#!/bin/bash
#$ -cwd

### workflow:
#  map reads using tophat (via bowtie)
#  estimate FPKM using cufflinks with reference genes

USAGE="Usage: $0 -i fastq -s setting [-n num_threads]"

nt=1 # default 1 threads for tophat
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


 bam=$outdir"/accepted_hits.bam"

## do cufflinks
cuffout=$bam"_cufflinks_ref"
# cuffout2=$bam"_cufflinks_ref-guide"

cmd="cufflinks -o $cuffout --GTF $GENES $bam"
# cmd="cufflinks -o $cuffout --compatible-hits-norm --GTF  $GENES $bam"
echo -e "do cufflinks with ref genes: \n $cmd"
$cmd

cuffout2=$bam"_cufflinks_ref-guide"
cmd2="cufflinks -o $cuffout2 --GTF-guide  $GENES $bam"
echo -e "do cufflinks with ref genes -guide: \n $cmd2"
$cmd2

qsub -l mem=2G,time=5::  /ifs/scratch/c2b2/ngs_lab/xs2182/code/getSNPs.sh $outdir $GENO

sh /ifs/scratch/c2b2/ngs_lab/xs2182/code/getCounts.sh $outdir"/accepted_hits.sam"  $outdir"/accepted_hits_counts.txt" $bam $GENES

qsub -l mem=1G,time=10:: /ifs/scratch/c2b2/ngs_lab/xs2182/code/QCplot.sh $outdir

# qsub -l mem=1G,time=5:: /ifs/scratch/c2b2/ngs_lab/xs2182/code/getCounts.sh $outdir"/accepted_hits.sam"  $outdir"/accepted_hits_counts.txt" $bam $GENES


