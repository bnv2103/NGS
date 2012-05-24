#!bin/bash
#$ -cwd
#$ -l mem=2G,time=4::


## New merge BAM script to minimize typo error, and validates existant of all files.

#  $1 Path/Project 
# $2 SAMPLE-NAME
# $3 $4 ... All RUN-ID to merge
	
SAMTOOLS="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/samtools"
BPATH="/ifs/data/c2b2/ngs_lab/ngs/code/NGS/Exome/"


proj=$1
shift
sample=$1
shift
numruns=$#
runs=$@
if [ $numruns -lt 2 ]
then
  echo "Usage: $0 PATH/project-id  SAMPLE-ID  RUN-ID1 RUN-ID2 ... RUN-IDN"
  echo "qsub -o log.sample-id.o -e log.sample-id.e mergeBAM.sh  /ifs/scratch/c2b2/ngs_lab/ngs/Projects/Exome-seq/120129_WENDY_WENDY_24_HUMAN_EXOME_60X_PE_HISEQ/ CDH632 120418_SN650_0270_AD0JGRACXX 120316_SN650_0258_AC0EAMACXX "
  exit
fi

cd $proj

#Sanity Check
flag=0
for arun in $runs
do
	if [[ ! -e "$arun/mapping/$sample.bam_refine/all.realigned.bam" ]]
	then
		echo "ERROR: $arun/mapping/$sample.bam_refine/all.realigned.bam does not exist ";
		flag=1
	fi
done
if [[ $flag == 1 ]]
then
	echo "ABORT: Some input file doesn't exist"
	exit
fi

#Create working directory

dirname="combined"
for arun in $runs
do
	x=`echo $arun | cut -f1 -d"_" `
	dirname=$dirname"_$x"
done
echo "working dir = $dirname "

if [ ! -e $dirname ]; then mkdir $dirname; fi
cd $dirname
if [ ! -e mapping ]; then mkdir mapping; fi
if [ ! -e logs ]; then mkdir logs;  fi
if [ ! -e global_setting.sh ]
then
	for arun in $runs
	do
		cat $proj/$arun/global_setting* >> global_setting.sh
	done
fi

setting="$proj/$dirname/global_setting.sh"
outfile="$proj/$dirname/mapping/$sample.bam"

sampledirname=`echo $outfile | sed 's/.bam//'`

count=0
suff=".bam_refine/all.realigned.bam"
infiles=""
for arun in $runs 
do
  count=`expr $count + 1`
  if [[ $count == 1 ]];then
      $SAMTOOLS view -H $proj/$arun/mapping/$sample$suff  > $sample.temp.header
	grep "@HD" $sample.temp.header > $sample.final.header
        grep "@SQ" $sample.temp.header >> $sample.final.header
	grep "@RG" $sample.temp.header  > $sample.rg.header
  else
  	$SAMTOOLS view -H $proj/$arun/mapping/$sample$suff |grep "@RG"  >> $sample.rg.header
  fi
  infiles=$infiles" $proj/$arun/mapping/$sample$suff "
done

echo "infiles = $infiles"
echo "settings = $setting "

cat $sample.rg.header >>$sample.final.header
grep "@PG" $sample.temp.header >> $sample.final.header
rm $sample.temp.header $sample.rg.header


$SAMTOOLS merge -f -h $sample.final.header $outfile $infiles
echo "Merge Complete"
$SAMTOOLS sort $outfile $outfile.sorted 
mv $outfile.sorted.bam $outfile
echo "Sort Complete"
$SAMTOOLS index $outfile
echo "Index Complete"
 

    OUTDIRNAME=$outfile"_refine"

    if [ ! -d $OUTDIRNAME ]; then
        mkdir $OUTDIRNAME
    fi
    OUTDIR=`readlink -f $OUTDIRNAME`
    status=$OUTDIR"/realign.status"
    if [ -e $status ] ; then
        rm -f $status
    fi
echo "Out directory = $OUTDIR "
echo "Output file = $outfile "
    touch $status
    mkdir -p $OUTDIR/logs/
    qmem=5 # default
    heapm=4
    qtime=24
    for (( i=1; i<=24; i++))
      do
      if [[ $i -lt 7 ]]; then  # mem=8 for large chr
          qmem=8
          heapm=7
          qtime=30
      fi

	if [[ $AUTO == "" ]];then
              cmd="qsub -N realign.$i.$sample -l mem=${qmem}G,time=${qtime}:: -o $OUTDIR/logs/realign.$i.o -e $OUTDIR/logs/realign.$i.e ${BPATH}/gatk_realign_atomic.sh -I $outfile -o $OUTDIR  -g $setting -L $i -c $status -m $heapm -A AUTO"
        else
              cmd="qsub -N realign.$i.$sample -l mem=${qmem}G,time=${qtime}:: -o $OUTDIR/logs/realign.$i.o -e $OUTDIR/logs/realign.$i.e ${BPATH}/gatk_realign_atomic.sh -I $output.bam -o $OUTDIR  -g $setting -L $i -c $status -m $heapm -A AUTO"
        fi
      echo $cmd
      $cmd
    done

