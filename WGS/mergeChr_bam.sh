#!/bin/bash
#$ -cwd
	
SAMTOOLS="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/samtools"
REF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/human_g1k_v37.fasta"
param=""
chr=$1
list=$2
outheader=$3
DIR=$4"/"

for dire in `cat $list` ;do
	param=$param$dire"/$chr "
done
CMD="$SAMTOOLS merge -f -h $DIR$outheader $DIR$chr.temp $param "
echo $CMD
$CMD
echo "Merge Complete"
$SAMTOOLS sort $DIR$chr.temp $DIR$chr.sorted
echo "Sort Complete on Merged BAM"
$SAMTOOLS index $DIR$chr.sorted.bam
echo "Index Complete on Merged BAM"
rm $DIR$chr.temp	

basedir=`dirname DIR`
LOG="$basedir/logs"
if [ ! -e $LOG ]
then
	mkdir $LOG
fi

CMD="qsub -o $LOG/$chr.refine.o -e $LOG/$chr.refine.e -N refine.$chr -l mem=10G,time=10::  /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/refine_bam.sh $chr $DIR "
echo $CMD
$CMD

exit

$SAMTOOLS rmdup $DIR$chr.sorted.bam $DIR$chr.sorted.bam.noDup.bam
echo "Remove PCR Duplicates complete"
# rm $DIR$chr.sorted.bam
## Should add -C int : Coefficient to cap mapping quality of poorly mapped reads. See the pileup command for details. [0]
$SAMTOOLS calmd -rEAb $DIR$chr.sorted.bam.noDup.bam $REF >  $DIR$chr.sorted.bam.noDup.bam.baq.bam
$SAMTOOLS index $DIR$chr.sorted.bam.noDup.bam.baq.bam
echo "Index Complete on bam.noDup.bam.baq.bam "
rm $DIR$chr.sorted.bam.noDup.bam
echo "Realigned using Samtools- CALMD with extended BAQ done"
# $SAMTOOLS sort -n  $DIR$chr.sorted.bam.noDup.bam.baq.bam $DIR$chr.sorted.bam.noDup.bam.baq.bam.namesorted
# echo "Name Sort Complete on bam.noDup.bam.baq.bam "
# $SAMTOOLS index $DIR$chr.sorted.bam.noDup.bam.baq.bam.sorted.bam
# echo "Index Complete on  bam.noDup.bam.baq.bam.sorted.bam "

 
qsub fixmte n stat
HEAP=7
TEMP=""
JAVA="java -Xmx${HEAP}g -Djava.io.tmpdir="${TEMP}
FIXMATE="$JAVA -jar ${PICARD}/FixMateInformation.jar"
CMD="$FIXMATE INPUT=$OUTDIR/$CHR.cleaned.bam OUTPUT=$OUTDIR/$CHR.fixed.bam SO=coordinate VALIDATION_STRINGENCY=SILENT "
echo "Fixmate: $CMD "
$CMD

if [[ $? == 0 ]]
    then
    echo "fixmate complete"
else
    echo "fixmate failed"
    exit 1
fi
$SAMTOOLS index $jsdfjahf.fixed.bam
echo "Index Complete on $jsdfjahf.fixed.bam "
# rm INPUT $jsdfjahf.fixed.bam

CMD="$JAVA -jar $FIXSTAT INPUT=${myroot}/bam/$a/refine/${i}.noDup.rl.sorted.fxmate.bam OUTPUT=${myroot}/bam/$a/refine/stat/${i}/${i}.noDup.rl.sorted.fxmate_details HISTOGRAM_FILE=${myroot}/bam/$a/refine/stat/${i}/${i}.noDup.rl.sorted.fxmate.hist.pdf REFERENCE_SEQUENCE=$REF "
echo "Fixstat: CollectInsertSizeMetrics \n $CMD "
$CMD

CMD="$JAVA -jar $GCbias INPUT=${myroot}/bam/$a/refine/${i}.noDup.rl.sorted.bam OUTPUT=${myroot}/bam/$a/refine/stat/${i}/${i}.noDup.rl.sorted.fxmate.gc.bias_detail CHART=${myroot}/bam/$a/refine/stat/${i}/${i}.gcBias.pdf REFERENCE_SEQUENCE=$REF "
echo "GC Bias: CollectGcBiasMetrics \n $CMD "
$CMD



