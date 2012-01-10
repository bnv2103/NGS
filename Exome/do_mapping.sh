#!/bin/bash
#$ -cwd

fq1list=$1
setting=$2

BPATH="/ifs/home/c2b2/ngs_lab/ngs/code/NGS/Exome"

if [[ $setting == "" ]]; then
    setting="/ifs/home/c2b2/ngs_lab/ngs/code/NGS/Exome/global_setting_b37.sh"
fi

. $setting 

mkdir -p mapping
mkdir -p logs
for f in `cat $fq1list`  ## list.forwardReads.fq.txt`
  do 
# example: chung32_lane_3_CDH-1-513_1.fastq
  fbase=`basename $f`	
  g=`echo $f | sed 's/1.fastq/3.fastq/'`
  sampleName=`basename $f | sed 's/\_1.fastq//' | cut -f4  -d '_'`
  
  cmd="qsub -pe smp 2  -l mem=4G,time=99:: -o logs/mapping.$fbase.o -e logs/mapping.$fbase.e -N map.$fbase $BPATH/mapping-two-cores.sh -i $f -p $g -z $sampleName -n $sampleName -s $setting -o mapping/$sampleName -t 2 -c 1"
  echo $cmd
  $cmd
done

