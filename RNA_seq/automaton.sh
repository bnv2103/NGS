#!/bin/bash
#$ -cwd


##  work flow: 
# 1. given a number, get the number of reads
# 2. run pipeline

USAGE="Usage: $0 -i fastq -s setting -n number-of-reads [-t num_threads] [-o outdir]"

nt=4 # default 4 threads for tophat
outdir=""

RUBY=`which ruby`

while getopts i:s:n:t:h opt
  do
  case "$opt" in
      i) fq="$OPTARG";;
      s) setting="$OPTARG";;
      n) num="$OPTARG";;
      t) nt="$OPTARG";;
      o) outdir="$OPTARG";;
      h) echo $USAGE; exit 1
  esac
done

if [[ $fq == "" || $setting == "" || $num == "" ]]
    then
    echo $USAGE
    exit 1
fi


## extract 


name=`basename $fq | sed 's/.fastq$//g' | sed 's/.fq$//g' | sed 's/\//\_/g'`
rfq=$name"_random_"$num.fq
cmd="$RUBY ~/code/NGS/utiles/randomly-take-N-reads-from-fastq.rb $fq $num > $rfq"

echo -e "randomly extract $num reads from $fq\n $cmd"
$cmd


## tophat and cufflinks

cmd="sh ~/code/NGS/RNA_seq/pipeline_cufflink-ref.sh  -i $rfq -s  $setting -n $nt "
echo -e "do tophat and cufflinks:\n $cmd"
$cmd
