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

## first do tophat

# get sample name
# g=`basename $fq | sed 's/.fastq//' | sed 's/\//\-/g'`
# rg=$g
# sampleName=$rg
# library=$rg


## SE
# cmd="tophat --phred64-quals --num-threads $nt  -o $outdir  --rg-id $rg --rg-sample $sampleName --rg-library $library $BOWTIEDB $fq"

# echo -e "do tophat: \n $cmd"
# $cmd

 bam=$outdir"/accepted_hits.bam"

# if [[ ! -s $bam ]]; then
#     echo "tophat failed"
#     exit 1
# fi


## do cufflinks
cuffout=$bam"_cufflinks_ref"
# cuffout2=$bam"_cufflinks_ref-guide"

cmd="cufflinks -o $cuffout --GTF $GENES $bam"
# cmd="cufflinks -o $cuffout --compatible-hits-norm --GTF  $GENES $bam"
echo -e "do cufflinks with ref genes: \n $cmd"
$cmd

## reference-guided assembly
# cmd2="cufflinks -o $cuffout2 --compatible-hits-norm --GTF-guide  $GENES $bam"
# cmd2="cufflinks -o $cuffout2 --GTF-guide  $GENES $bam"
# echo -e "do cufflinks with ref genes -guide: \n $cmd2"
# $cmd2


# for f in *cufflinks; do ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/comb_stats.rb $f "summary.csv"; done
# ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/comb_stats.rb $outdir $GENO"summary.csv"

if [[ $GENO == "mouse" ]];
    then
    qsub -l mem=2G,time=10:: /ifs/scratch/c2b2/ngs_lab/xs2182/code/getBed_mouse.sh $bam $outdir
    qsub -l mem=2G,time=5:: /ifs/scratch/c2b2/ngs_lab/xs2182/code/getSNPs.sh $outdir "mouse"
fi
if [[ $GENO == "human" ]];
    then
    qsub -l mem=2G,time=10:: /ifs/scratch/c2b2/ngs_lab/xs2182/code/getBed.sh $bam $outdir
    qsub -l mem=2G,time=5:: /ifs/scratch/c2b2/ngs_lab/xs2182/code/getSNPs.sh $outdir "human"
fi
if [[ $GENO == "rat" ]];
    then
    qsub -l mem=2G,time=10:: /ifs/scratch/c2b2/ngs_lab/xs2182/code/getBed_mouse.sh $bam $outdir
    qsub -l mem=2G,time=5:: /ifs/scratch/c2b2/ngs_lab/xs2182/code/getSNPs.sh $outdir "rat"
fi



# sh /ifs/scratch/c2b2/ngs_lab/xs2182/code/getSNPs.sh
# ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/mergeVPC.rb

