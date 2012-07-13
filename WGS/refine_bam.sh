#!/bin/bash
#$ -cwd
	
SAMTOOLS="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/samtools"
REF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/human_g1k_v37.fasta"
chr=$1
DIR=$2"/"


$SAMTOOLS rmdup $DIR$chr.sorted.bam $DIR$chr.sorted.bam.noDup.bam
echo "Remove PCR Duplicates complete"
# rm $DIR$chr.sorted.bam
## Should add -C int : Coefficient to cap mapping quality of poorly mapped reads. See the pileup command for details. [0]
$SAMTOOLS calmd -rEAb $DIR$chr.sorted.bam.noDup.bam $REF >  $DIR$chr.sorted.bam.noDup.bam.baq.bam
$SAMTOOLS index $DIR$chr.sorted.bam.noDup.bam.baq.bam
echo "Index Complete on bam.noDup.bam.baq.bam "
rm $DIR$chr.sorted.bam.noDup.bam
echo "Realigned using Samtools- CALMD with extended BAQ done"

exit

TEMP=$DIR$chr"_temp"
if [ ! -e $TEMP ];then
	mkdir $TEMP
fi

JAVA="java -Xmx7g -Djava.io.tmpdir=$TEMP"
FIXMATE="$JAVA -jar ${PICARD}/FixMateInformation.jar"
$FIXMATE INPUT=$OUTDIR/$CHR.cleaned.bam OUTPUT=$OUTDIR/$CHR.fixed.bam SO=coordinate VALIDATION_STRINGENCY=SILENT


out=${DIR}/$chr.fixing.sh
        echo '#!/bin/bash'  > $out
        echo 'uname -a' >> $out
        echo "source $myroot/$setting" >> $out
        echo "java -Xmx2g -Djava.io.tmpdir=${dir}/temp -jar $FIXMATE INPUT=${dir}/bam/$a/refine/$bam OUTPUT=${myroot}/bam/$a/refine/${i}.noDup.rl.sorted.fxmate.bam SO=coordinate VALIDATION_STRINGENCY=SILENT" >> $out
        cmd="java -Xmx2g -Djava.io.tmpdir=${dir}/temp -jar $FIXSTAT INPUT=${dir}/bam/$a/refine/${i}.noDup.rl.sorted.fxmate.bam OUTPUT=${myroot}/bam/$a/refine/stat/${i}/${i}.noDup.rl.sorted.fxmate_details HISTOGRAM_FILE=${myroot}/bam/$a/refine/stat/${i}/${i}.noDup.rl.sorted.fxmate.hist.pdf REFERENCE_SEQUENCE=$REF"
        echo $cmd >> $out
        cmd="java -Xmx2g -Djava.io.tmpdir=${dir}/temp -jar $GCbias INPUT=${dir}/bam/$a/refine/${i}.noDup.rl.sorted.bam OUTPUT=${myroot}/bam/$a/refine/stat/${i}/${i}.noDup.rl.sorted.fxmate.gc.bias_detail CHART=${myroot}/bam/$a/refine/stat/${i}/${i}.gcBias.pdf REFERENCE_SEQUENCE=$REF"
        echo $cmd >> $out
        cmd="samtools index ${dir}/bam/$a/refine/${i}.noDup.rl.sorted.fxmate.bam"
        echo $cmd >> $out
        echo 'kill -USR2 $watch_pid; wait' >> $out
        qsub -l mem=3G,time=2:: -o ${dir}/$log/fm-${i}.o -e ${dir}/$log/fm-${i}.e $out
 

