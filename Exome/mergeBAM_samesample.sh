#!bin/bash
#$ -cwd
#$ -l mem=2G,time=4::

##give full path
	
SAMTOOLS="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/samtools"
BPATH="/ifs/data/c2b2/ngs_lab/ngs/code/NGS/Exome/"
setting=$1
shift
outfile=$1
shift
if [ $# -lt 2 ]
then
  echo "Usage: $0 PATH/Global_setting.sh PATH/OUTPUT.bam  inp1.realigned.bam inp2.realigned.bam  .. inpN"
  echo "qsub -o log.cdh632.o -e log.cdh632.e mergeBAM.sh  /ifs/scratch/c2b2/ngs_lab/ngs/Projects/Exome-seq/120129_WENDY_WENDY_24_HUMAN_EXOME_60X_PE_HISEQ/combineReads120418_120316/global_setting.sh   /ifs/scratch/c2b2/ngs_lab/ngs/Projects/Exome-seq/120129_WENDY_WENDY_24_HUMAN_EXOME_60X_PE_HISEQ/combineReads120418_120316/mapping/CDH632.bam /ifs/scratch/c2b2/ngs_lab/ngs/Projects/Exome-seq/120129_WENDY_WENDY_24_HUMAN_EXOME_60X_PE_HISEQ/120418_SN650_0270_AD0JGRACXX/mapping/CDH632.bam_refine/all.realigned.bam /ifs/scratch/c2b2/ngs_lab/ngs/Projects/Exome-seq/120129_WENDY_WENDY_24_HUMAN_EXOME_60X_PE_HISEQ/120316_SN650_0258_AC0EAMACXX/mapping/CDH632.bam_refine/all.realigned.bam "
  exit
fi

dirname=`echo $outfile | sed 's/.bam//'`
if [[ ! -e $dirname ]]; then mkdir $dirname; fi
cd $dirname

nofiles=$#
infiles=$@
if [ -e temp.header ]; then rm temp.header ; fi
if [ -e rg.header ]; then rm rg.header ; fi

count=0
# infiles_final=""
for i in $infiles 
do
  count=`expr $count + 1`
  idir=`dirname $i `
  if [[ $count == 1 ]];then
      $SAMTOOLS view -H  $i > temp.header
	grep "@HD" temp.header > final.header
        grep "@SQ" temp.header >> final.header
	grep "@RG" temp.header  > rg.header
  else
  	$SAMTOOLS view -H  $i |grep "@RG"  >> rg.header
  fi
#  mv $i $idir"/"$count".bam"
#  infiles_final="$infiles_final $idir/$count.bam"
done
	cat rg.header >>final.header
	grep "@PG" temp.header >> final.header
rm temp.header rg.header

$SAMTOOLS merge -f -h final.header $outfile $infiles
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
glodir=`dirname $OUTDIR`        
glodir1=`dirname $glodir`
echo $OUTDIR
echo $setting
echo $outfile
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

      g=`basename $outfile | sed 's/\//_/g'`
if [[ $AUTO == "" ]];then
              cmd="qsub -N realign.$i.$g -l mem=${qmem}G,time=${qtime}:: -o $OUTDIR/logs/realign.$i.o -e $OUTDIR/logs/realign.$i.e ${BPATH}/gatk_realign_atomic.sh -I $outfile -o $OUTDIR  -g $setting -L $i -c $status -m $heapm -A AUTO"
        else
              cmd="qsub -N realign.$i.$g -l mem=${qmem}G,time=${qtime}:: -o $OUTDIR/logs/realign.$i.o -e $OUTDIR/logs/realign.$i.e ${BPATH}/gatk_realign_atomic.sh -I $output.bam -o $OUTDIR  -g $setting -L $i -c $status -m $heapm -A AUTO"
        fi
      echo $cmd
      $cmd
    done

