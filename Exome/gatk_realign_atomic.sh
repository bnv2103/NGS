#!/bin/bash
#$ -cwd
# Findmem

# Realign reads locally around indels
# Run with all samples from a single chromosome together in $INP.list

HEAP=4
### TEMP="/ifs/scratch/c2b2/af_lab/saec/temp/"
INP=""
CHR=""
MEM=""
chain=""
OUTDIR=""
MaxReads=20000
USAGE="Usage: $0 -I <Input bam file> -L <Chromosome> -g <global config> [-m heap_memory]"
ERRORMESSAGE="#### ERROR"
ERRORMESSAGE1="The following error has occurred"

while getopts I:L:m:g:r:o:c:h o
  do 
  case "$o" in
      I) INP="$OPTARG";;
      L) CHR="$OPTARG";;
      m) MEM="$OPTARG";;
      o) OUTDIR="$OPTARG";;
      r) MaxReads="$OPTARG";;
      c) chain="$OPTARG";;  # if $chain != 0,  chain reaction.
      g) GLOBAL="$OPTARG";;  # global config
      h) echo $USAGE
	  exit 1;;
  esac
done


if [[ $INP == ""  || $GLOBAL == "" ]]
then
        echo $USAGE
        exit 1
fi

. $GLOBAL

if [ ! $MEM == "" ]
    then
    HEAP=$MEM
fi

if [[ $OUTDIR == "" ]]; then
    OUTDIR=$INP"_refine"
fi

if [ ! -d $OUTDIR ]; then
    mkdir $OUTDIR
fi

if [ ! $chain == "" ]; then
    status=$chain
else
    status=$OUTDIR"/realign.status"
fi

touch $status


if [ $JOB_ID == "" ]; then
    JOB_ID=$CHR
fi
 
# JOB_ID is the qsub job ID
TEMP=$OUTDIR"/temp/"$JOB_ID"/"

if [ ! -d $TEMP ]; then
    mkdir -p $TEMP
fi

JAVA="java -Xmx${HEAP}g -Djava.io.tmpdir="${TEMP}
GATK="$JAVA -jar "${GATKJAR}
FIXMATE="$JAVA -jar ${PICARD}/FixMateInformation.jar" 

if [[ $CHR == "" ]]
    then
    echo "no chr prodived"
    exit 1
fi
if [[ $CHR == "23" ]]
    then
    CHR="X"
fi
if [[ $CHR == "24"  ]]
    then
    CHR="Y"
fi
if [[ $REFTYPE == "hg" ]]  # hg18/19 -> chr1, chr2 etc;  build36/37 -> 1, 2 etc
    then
    CHR="chr${CHR}"
fi

echo $CHR 

$GATK \
    -L $CHR \
    -T RealignerTargetCreator \
    -I $INP \
    -R $REF \
    -D $DBSNP \
    -DBQ 10 \
    -B:indels,VCF $INDELVCF  \
    -o $OUTDIR/$CHR.forRealigner.intervals

$GATK \
    -L $CHR \
    -I $INP \
    -R $REF \
    -D $DBSNP \
    -DBQ 10 \
    -T IndelRealigner \
    -B:indels,VCF $INDELVCF  \
    --maxReadsForRealignment $MaxReads \
    -compress 5 \
    -targetIntervals $OUTDIR/$CHR.forRealigner.intervals \
    -o $OUTDIR/$CHR.cleaned.bam 

# Need to implement for known-only/lane level cleaning?
if [[ $? == 0 ]]
    then
    echo "local realign complete"
else
    echo "local realign failed"
    exit 1
fi

$FIXMATE INPUT=$OUTDIR/$CHR.cleaned.bam OUTPUT=$OUTDIR/$CHR.fixed.bam SO=coordinate VALIDATION_STRINGENCY=SILENT 


if [[ $? == 0 ]]
    then
    rm -f $OUTDIR/$CHR.cleaned.bam*
    echo "fixmate complete"
    echo $CHR >> $status
else
    echo "fixmate failed"
    exit 1
fi
# rm -rf $OUTDIR/$CHR.realign.*



### "chain reaction": check the status of other chromosomes. If all chromosomes are complete, then proceed to the next step (merge and recalibrate etc). 
if [[ $chain != "" ]]; then
# wait for 6 seconds if not complete
    for (( i=0; i<2; i++ ))
      do 
      completed=`wc -l $status | awk '{print $1}'`
      if [[ $completed == "24" ]]; then  # all chromosomes completed, proceed to the next step (merge), 
	  echo "all completed" >>  $status  # make it stop
	  
	  if [ ! -e $OUTDIR"/forRealigner.intervals.bz2" ]
	      then
	      cat $OUTDIR/*.forRealigner.intervals | bzip2 - > $OUTDIR/forRealigner.intervals.bz2
	      rm -f $OUTDIR/*.forRealigner.intervals
	  fi
	  
##  merge bam files
	  if [ ! -e $OUTDIR"/all.realigned.bam"  ] # need merge
	      then
	      echo "merge bam files"
	      ${SAMTOOLS} merge $OUTDIR/all.realigned.bam $OUTDIR/*.fixed.bam  > $OUTDIR/log.merge 2>&1
	      
	      rm -f $OUTDIR/*.fixed.bam 
	      ${SAMTOOLS} index $OUTDIR/all.realigned.bam  > $OUTDIR/log.index 2>&1

	      if [[  -s $OUTDIR/all.realigned.bam  ]]; then
		  rm -rf $INP  # delete original bam

		  # trigger downstream analysis (dedup and recalibrate)
		  echo "dedup"
		  
		  cmd="${BPATH}/deduplication.sh -i $OUTDIR/all.realigned.bam -s $GLOBAL"
		  echo $cmd
		  $cmd
		  
		  echo "recalibrate quality"
		  date 

		  mkdir -p $OUTDIR/logs/
		  cmd="qsub -l mem=8G,time=55:: -o $OUTDIR/logs/all.realigned.bam.log.recalib.o -e $OUTDIR/logs/all.realigned.bam.log.recalib.e  ${BPATH}/gatk_recalibrate.sh -m 6000  -g $GLOBAL   -I $OUTDIR/all.realigned.bam  -o $OUTDIR/all.recalibrated.bam -c 1"
		  echo $cmd
		  $cmd
	      fi
	      


	  fi
	  
	  i=2
	  
      else
	  sleep 3
      fi
    done
fi


rm -fr $TEMP
