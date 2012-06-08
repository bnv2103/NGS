#!/bin/bash
#$ -cwd

fqdir=$1	#Takes single directory name of the sample filename not a lit	
setting=$2
automated=$3	#if initiated by automated pipeline then input argument must be ProjectID, this triggers automated downstream steps

BPATH="/ifs/data/c2b2/ngs_lab/ngs/code/NGS/Exome"
BPATH="/ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/"
if [[ $setting == "" ]]; then
    setting="/ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/global_setting_b37.sh"
fi

. $setting 

if [[ !  -d mapping ]] ;then  mkdir mapping; fi
if [[ !  -d logs ]] ;then  mkdir logs; fi

# example: R_U_N_ID_laneX_SAMPLEID_1.fastq   = 120210_SN650_0253_AD0JGWACXX_lane4_Misc-1418_1.fastq
for fq in `ls $fqdir/x* `; do
  job_name=`basename $fqdir |  sed 's/\_1.fastq//' `	#unique name for jobname (includes runid_lane#_sampleID
  job_suff=`basename $fq`
  g=`echo $fq | sed 's/_1.fastq/_3.fastq/'`
  nthread=2
  if  [ ! -e $g ];  then g="" ;  nthread=2 ;else   nthread=4; fi
  sampleName=` echo $job_name | cut -f6  -d '_'`	#acc to  new convention, to change to old do -f4
  lane=` echo $job_name | cut -f5  -d '_'`
  job_id=$job_name"_"$lane"_"$job_suff

  if [[ ! -d mapping/$sampleName ]] ;then  mkdir mapping/$sampleName; fi
  if [[ ! -d mapping/$sampleName/$lane ]] ;then  mkdir mapping/$sampleName/$lane; fi

   rg1=`basename $fqdir|cut -f3 -d"_"`
   ln=`echo $lane | sed 's/lane//'`
   readgroup=$rg1"_"$ln"_"$sampleName
   ID="$rg1$ln$sampleName"

if [[ $automated == "" ]]; #was NOT triggered by automatic pipeline
then
  cmd="qsub -pe smp $nthread -R y -l mem=4G,time=48:: -o logs/mapping.$job_id.o -e logs/mapping.$job_id.e -N map.$job_id $BPATH/mapping-two-cores.sh -i $fq -p $g -y $ID -z $readgroup -n $sampleName -s $setting -o mapping/$sampleName/$lane/$job_suff -t 2 -c 1 "
else	#trigger automatic process
  cmd="qsub -pe smp $nthread -R y -l mem=4G,time=48:: -o logs/mapping.$job_id.o -e logs/mapping.$job_id.e -N map.$job_id.AUTO $BPATH/mapping-two-cores.sh -i $fq -p $g -y $ID  -z $readgroup -n $sampleName -s $setting -o mapping/$sampleName/$lane/$job_suff -t 2 -c 1 -A $automated"
fi
  echo $cmd >> logs/history.$sampleName.$lane.txt
  $cmd
done

