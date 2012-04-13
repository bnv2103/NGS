#!/bin/bash
#$ -cwd

fq=$1	#Takes single filename not a lit	
setting=$2
automated=$3	#if initiated by automated pipeline then input argument must be ProjectID, this triggers automated downstream steps

BPATH="/ifs/data/c2b2/ngs_lab/ngs/code/NGS/Exome"

if [[ $setting == "" ]]; then
    setting="/ifs/home/c2b2/ngs_lab/ngs/code/NGS/Exome/global_setting_b37.sh"
fi

. $setting 

#for f in `cat $fq1list`  ## list.forwardReads.fq.txt`
#  do 
# example: RUNID_laneX_SAMPLEID_1.fastq
  job_id=`basename $fq |  sed 's/\_1.fastq//' `	#unique name for jobname (includes runid_lane#_sampleID
  g=`echo $fq | sed 's/_1.fastq/_3.fastq/'`
  nthread=2
  if  [ ! -e $g ];  then g="" ;  nthread=2 ;else   nthread=4;fi
  sampleName=`basename $fq | sed 's/\_1.fastq//' | cut -f6  -d '_'`	#acc to  new convention, to change to old do -f4
  
if [[ $automated == "" ]]; #was NOT triggered by automatic pipeline
then
  cmd="qsub -pe smp $nthread -R y -l mem=4G,time=48:: -o logs/mapping.$sampleName.o -e logs/mapping.$sampleName.e -N map.$job_id $BPATH/mapping-two-cores.sh -i $fq -p $g -z $sampleName -n $sampleName -s $setting -o mapping/$sampleName -t 2 -c 1 "
else	#trigger automatic process
  cmd="qsub -pe smp $nthread -R y -l mem=4G,time=48:: -o logs/mapping.$sampleName.o -e logs/mapping.$sampleName.e -N map.$job_id.AUTO $BPATH/mapping-two-cores.sh -i $fq -p $g -z $sampleName -n $sampleName -s $setting -o mapping/$sampleName -t 2 -c 1 -A $automated"
fi
  echo $cmd >> logs/history.$sampleName.txt
  $cmd
#done
