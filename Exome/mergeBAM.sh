#!bin/bash
#$ -cwd
#$ -l mem=2G,time=4::

##give full path
	
SAMTOOLS="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/samtools"
BPATH="/ifs/data/c2b2/ngs_lab/ngs/code/NGS/Exome/"
setting=$1
shift
outfile=$2
shift
if [ $# -lt 2 ]
then
  echo "Usage: $0 arg1 arg2 .. argN"
  exit
fi

dirname=`echo $outfile | sed 's/.bam//'`
if [[ ! -e $dirname ]]; then mkdir $dirname; fi
cd $dirname

nofiles=$#
infiles=$@
if [ -e all.header ]; then rm all.header ; fi
if [ -e rg.header ]; then rm rg.header ; fi
count=0
infiles_final=""
for i in $infiles 
do
  count=`expr $count + 1`
#  if [ -e $count.bam ]; then rm $count.bam; fi
#  ln -s $i $count.bam
  idir=`dirname $i `
  if [[ $count == 1 ]];then
      $SAMTOOLS view -H  $i > temp.header
	grep "@HD" temp.header > final.header
        grep "@SQ" temp.header >> final.header
	grep "@RG" temp.header |sed "s/ID:.\+\tPL/ID:$count\tPL/" >> rg.header
# grep "@RG" temp.header |sed "s/ID:/$count/
  else
  	$SAMTOOLS view -H  $i |grep "@RG"|sed "s/ID:.\+\tPL/ID:$count\tPL/"  >> rg.header
  fi
  mv $i $idir"/"$count".bam"
  infiles_final.=" $idir/$count.bam"
done
	cat rg.header >>final.header
	grep "@PG" temp.header >> final.header
rm temp.header

$SAMTOOLS merge -f -r -h final.header $outfile $infiles_final
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
#      setting=`ls $glodir1"/global_setting*" `
echo $OUTDIR
echo $setting
echo $outfile
# exit
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
              cmd="qsub -N realign.$i.$g -l mem=${qmem}G,time=${qtime}:: -o $OUTDIR/logs/realign.$i.o -e $OUTDIR/logs/realign.$i.e ${BPATH}/gatk_realign_atomic.sh -I $outfile -o $OUTDIR  -g $setting -L $i -c $status -m $heapm "
        else
              cmd="qsub -N realign.$i.$g -l mem=${qmem}G,time=${qtime}:: -o $OUTDIR/logs/realign.$i.o -e $OUTDIR/logs/realign.$i.e ${BPATH}/gatk_realign_atomic.sh -I $output.bam -o $OUTDIR  -g $setting -L $i -c $status -m $heapm -A AUTO"
        fi
      echo $cmd
      $cmd
    done

