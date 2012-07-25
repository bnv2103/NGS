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

# cmd="tophat --phred64-quals --num-threads $nt  -o $outdir  --rg-id $rg --rg-sample $sampleName --rg-library $library --mate-inner-dist 100 $BOWTIEDB $fq1 $fq2"

# echo -e "do tophat: \n $cmd"
# $cmd

# tophat -o tophat_MCF7 -p 8 --fusion-search --keep-fasta-order --bowtie1 --no-coverage-search -r 0 --mate-std-dev 80 --fusion-min-dist 100000 --fusion-anchor-length 13 --fusion-ignore-chromosomes chrM /path/to/h_sapiens/bowtie_index SRR064286_1.fastq SRR064286_2.fastq 

cmd="tophat --phred64-quals -o $outdir --num-threads $nt --rg-id $rg --rg-sample $sampleName --rg-library $library --keep-fasta-order --no-coverage-search -r 50 --mate-std-dev 100  --fusion-min-dist 100000 --fusion-anchor-length 13 --fusion-ignore-chromosomes chrM $BOWTIEDB $fq1 $fq2"

echo -e "do tophat: \n $cmd"
$cmd
 
# tophat-fusion-post -p 8 --num-fusion-reads 1 --num-fusion-pairs 2 --num-fusion-both 5 /path/to/h_sapiens/bowtie_index

cmd="tophat-fusion-post -p $nt --num-fusion-reads 1 --num-fusion-pairs 2 --num-fusion-both 5 $BOWTIEDB" 

bam=$outdir"/accepted_hits.bam"

if [[ ! -s $bam ]]; then
    echo "tophat failed"
    exit 1
fi


## do cufflinks
cuffout=$bam"_cufflinks_ref"

# cmd="cufflinks -o $cuffout --compatible-hits-norm --GTF  $GENES $bam"
cmd="cufflinks -o $cuffout --GTF $GENES $bam"
echo -e "do cufflinks with ref genes: \n $cmd"
$cmd

## do cufflinks with guide
cuffout=$bam"_cufflinks_ref-guide"

cmd1="cufflinks -o $cuffout --GTF-guide  $GENES $bam"

echo -e "do cufflinks without ref genes: \n $cmd1"
$cmd1

# after run cufflink, statistics of number of reads and FPKM
ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/comb_stats_PE.rb $outdir "summary.csv"

cd $outdir
ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/QCplot.rb accepted_hits.bam "QC.pdf"
cd ..

samtools sort $bam $outdir"/accepted_hits.sorted"
samtools index $outdir"/accepted_hits.sorted.bam"



 if [[ $GENO == "mouse" ]];
    then
    qsub -l mem=2G,time=5:: -o logs/getSNPs.o -e logs/getSNPs.e /ifs/scratch/c2b2/ngs_lab/xs2182/code/getSNPs.sh $outdir "mouse"
 fi
 if [[ $GENO == "human" ]];
    then
    qsub -l mem=2G,time=5:: -o logs/getSNPs.o -e logs/getSNPs.e /ifs/scratch/c2b2/ngs_lab/xs2182/code/getSNPs.sh $outdir "human"
 fi

