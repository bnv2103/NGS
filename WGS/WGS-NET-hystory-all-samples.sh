##remove PCR duplicate
setting='global_setting_b37.sh'
 . $setting
Inputdir="/ifs/scratch/c2b2/ngs_lab/ngs/Projects/WG-seq/120309_ANDREA_GABRIELLE_6_HUMAN_GENOME_30X_PE_HISEQ/SAMPLES"
sampleList="sample.list.all.2"
bamlist="bam.name.list"
mkdir -p bam
dir=`pwd`
#from bam file elimina PCR duplicate read with samtools for each sample for each chromosome
for a in `cat $sampleList`
do
	mkdir -p bam/$a
	mkdir -p bam/$a/log
	log="bam/$a/log"
	for bam in `cat $bamlist`
	do
	echo $bam
	i=`echo $bam $|cut -d "." -f1 `
	echo $i
	file="${dir}/bam/$a/${i}.noDup.bam"
	num=`grep finish ${dir}/$log/removeDupl-${i}}.o | wc -l`
	echo $file
	if [[ -e "$file" && "$num" -eq 1 ]]
	then
	echo "here $file"
	else
	rm ${dir}/bam/$a/log/removeDupl-${i}*
	out=${dir}/bam/$a/log/removeDupl-${i}.sh
	echo '#!/bin/bash'  > $out
	echo 'uname -a' >> $out
	echo "source $myroot/$setting" >> $out
	echo 'watch-myself &' >> $out
	cmd="samtools rmdup  ${Inputdir}/$a/BAM/$bam ${dir}/bam/$a/${i}.noDup.bam "
	echo $cmd >> $out
	echo 'kill -USR2 %1' >> $out
	qsub -l mem=4G,time=5:: -o ${dir}/$log/removeDupl-${i}.o -e ${dir}/$log/removeDupl-${i}.e $out
	fi
	done
done

##called is the samtool to reallined
setting='global_setting_b37.sh'
 . $setting
sampleList="sample.list.all.2"
bamlist="bam.name.list"
dir=`pwd`
for a in `cat $sampleList`
do
	mkdir -p bam/$a
	mkdir -p bam/$a/refine
	mkdir -p bam/$a/refine/log
	log="bam/$a/refine/log"
	/bin/ls -f bam/$a/*bam > $a.bam.noDup.list
	for al in `cat $a.bam.noDup.list`
	do
	bam=`basename $al`
	i=`echo $bam |cut -d "." -f1 `
	echo $i
	file="${dir}/bam/$a/refine/${i}.noDup.rl.bam"
	num=`grep finish ${dir}/$log/calmd-${i}.o | wc -l`
	echo $file
	if [[ -e "$file" && "$num" -eq 1 ]]
	then
	echo "here $file"
	else
	rm ${dir}/bam/$a/refine/log/calmd-${i}* ##attention it is a rm option
	out=${dir}/$log/calmd-${i}.sh
	echo '#!/bin/bash'  > $out
	echo 'uname -a' >> $out
	echo 'watch-myself &' >> $out
	echo "source $myroot/$setting" >> $out
	cmd="samtools calmd -rEAb ${dir}/bam/${a}/${bam} $REF > ${dir}/bam/$a/refine/${i}.noDup.rl.bam"
	echo $cmd >> $out
	echo 'kill -USR2 %1' >> $out
	echo $cmd
	qsub -l mem=4G,time=10:: -o ${dir}/$log/calmd-${i}.o -e ${dir}/$log/calmd-${i}.e $out
	fi
	done
done


##index and sorting
#
setting='global_setting_b37.sh'
 . $setting
sampleList="sample.list.all.all"
dir=`pwd`
for a in `cat $sampleList`
do
	#/bin/ls bam/$a/refine/*.noDup.rl.bam > $a.noDup.rl.bam.list
	for al in `cat $a.noDup.rl.bam.list` 
	do
	bam=`basename $al`
	i=`echo $bam |cut -d "." -f1 `
	#echo $i
	mkdir -p bam/$a/refine/stat
	mkdir -p bam/$a/refine/stat/${i}
	log="bam/$a/refine/log"
	num=`grep finish ${dir}/$log/ixSrSt-${i}.o | wc -l`
	file="${dir}/bam/$a/refine/$i.noDup.rl.sorted.bam"
	if [[ -e $file && "$num" -eq 1 ]]
	then
	continue #echo "here $file $num"
	else
	echo "$file $num"
	if [[ "$i" =~ "GL000" ]]
	then
	echo "matched"
	#rm ${dir}/$log/ixSrSt-${i}.*
	out=${dir}/$log/ixSrSt-${i}.sh
	echo '#!/bin/bash'  > $out
	echo 'uname -a' >> $out
	echo 'watch-myself &' >> $out
	echo "source $myroot/$setting" >> $out
	cmd="samtools sort ${dir}/bam/$a/refine/$bam ${dir}/bam/$a/refine/$i.noDup.rl.sorted"
	echo $cmd >> $out
	cmd="samtools index ${dir}/bam/$a/refine/$i.noDup.rl.sorted.bam"
	echo $cmd >> $out
	cmd="samtools idxstats ${dir}/bam/$a/refine/$i.noDup.rl.sorted.bam > ${dir}/bam/$a/refine/stat/${i}/$i.$a.noDup.rl.idxstats"
	echo $cmd >> $out
	cmd="samtools flagstat ${dir}/bam/$a/refine/$i.noDup.rl.sorted.bam > ${dir}/bam/$a/refine/stat/${i}/${i}.${a}.noDup.rl.flagstat"
	echo $cmd >> $out
	echo "rm ${dir}/bam/$a/${i}.noDup.bam " >> $out
	echo 'kill -USR2 %1' >> $out
	echo $out
	#qsub -l mem=2G,time=2:: -o ${dir}/$log/ixSrSt-${i}.o -e ${dir}/$log/ixSrSt-${i}.e $out
	else	
	#rm ${dir}/$log/ixSrSt-${i}*
	out=${dir}/$log/ixSrSt-${i}.sh
	echo '#!/bin/bash'  > $out
	echo 'uname -a' >> $out
	echo 'watch-myself &' >> $out
	echo "source $myroot/$setting" >> $out
	cmd="samtools sort ${dir}/bam/$a/refine/$bam ${dir}/bam/$a/refine/$i.noDup.rl.sorted"
	echo $cmd >> $out
	cmd="samtools index ${dir}/bam/$a/refine/$i.noDup.rl.sorted.bam"
	echo $cmd >> $out
	cmd="samtools idxstats ${dir}/bam/$a/refine/$i.noDup.rl.sorted.bam > ${dir}/bam/$a/refine/stat/${i}/$i.$a.noDup.rl.idxstats"
	echo $cmd >> $out
	cmd="samtools flagstat ${dir}/bam/$a/refine/$i.noDup.rl.sorted.bam > ${dir}/bam/$a/refine/stat/${i}/${i}.${a}.noDup.rl.flagstat"
	echo $cmd >> $out
	echo "rm ${dir}/bam/$a/${i}.noDup.bam " >> $out
	echo 'kill -USR2 %1' >> $out
	echo $out
	qsub -l mem=8G,time=8:: -o ${dir}/$log/ixSrSt-${i}.o -e ${dir}/$log/ixSrSt-${i}.e $out
	fi
	fi
	done
done 

##fixmate and stat on mate
setting='global_setting_b37.sh'
 . $setting
sampleList="sample.list.all.all"
dir=`pwd`
for a in `cat $sampleList`
do
	mkdir -p bam/$a
	mkdir -p bam/$a/refine
	mkdir -p bam/$a/refine/log
	mkdir -p bam/$a/refine/stat/
	log="bam/$a/refine/log"
	mkdir -p temp
	/bin/ls bam/$a/refine/*.noDup.rl.sorted.bam > $a.noDup.rl.bam.sorted.list
	for al in `cat $a.noDup.rl.bam.sorted.list`
	do
	bam=`basename $al`
	i=`echo $bam |cut -d "." -f1 `
	echo $i
	mkdir -p bam/$a/refine/stat/${i}
	file="${dir}/bam/$a/refine/stat/${i}/${i}.noDup.rl.sorted.fxmate_details"
	echo $file
	num=`grep finish ${dir}/$log/fm-${i}.o | wc -l`
	if [[ -e $file && "$num" -eq 1 ]]
	then
	echo "here $file $num"
	else
	echo "pippo $num"
	if [[ "$i" =~ "GL000" ]]
	then
	echo "matched"
	rm ${dir}/$log/fm-${i}.*
	out=${dir}/$log/fm-${i}.sh
	echo '#!/bin/bash'  > $out
	echo 'uname -a' >> $out
	echo 'watch-myself & watch_pid=$!' >> $out
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
	else
	rm ${dir}/$log/fm-${i}.*
	out=${dir}/$log/fm-${i}.sh
	echo '#!/bin/bash'  > $out
	echo 'uname -a' >> $out
	echo 'watch-myself & watch_pid=$!' >> $out
	echo "source $myroot/$setting" >> $out
	echo "java -Xmx7g -Djava.io.tmpdir=${dir}/temp -jar $FIXMATE INPUT=${dir}/bam/$a/refine/$bam OUTPUT=${myroot}/bam/$a/refine/${i}.noDup.rl.sorted.fxmate.bam SO=coordinate VALIDATION_STRINGENCY=SILENT" >> $out
	cmd="java -Xmx7g -Djava.io.tmpdir=${dir}/temp -jar $FIXSTAT INPUT=${dir}/bam/$a/refine/${i}.noDup.rl.sorted.fxmate.bam OUTPUT=${myroot}/bam/$a/refine/stat/${i}/${i}.noDup.rl.sorted.fxmate_details HISTOGRAM_FILE=${myroot}/bam/$a/refine/stat/${i}/${i}.noDup.rl.sorted.fxmate.hist.pdf REFERENCE_SEQUENCE=$REF"
	echo $cmd >> $out
	cmd="java -Xmx7g -Djava.io.tmpdir=${dir}/temp -jar $GCbias INPUT=${dir}/bam/$a/refine/${i}.noDup.rl.sorted.bam OUTPUT=${myroot}/bam/$a/refine/stat/${i}/${i}.noDup.rl.sorted.fxmate.gc.bias_detail CHART=${myroot}/bam/$a/refine/stat/${i}/${i}.gcBias.pdf REFERENCE_SEQUENCE=$REF"
	echo $cmd >> $out
	cmd="samtools index ${dir}/bam/$a/refine/${i}.noDup.rl.sorted.fxmate.bam"
	echo $cmd >> $out
	echo 'kill -USR2 $watch_pid; wait' >> $out
	qsub -l mem=8G,time=8:: -o ${dir}/$log/fm-${i}.o -e ${dir}/$log/fm-${i}.e $out 
	fi
	fi
	done
done

####################
###qc after reallining GARK
###################

setting='global_setting_b37.sh'
 . $setting
sampleList="sample.list.all.2"
dir=`pwd`
for a in `cat $sampleList`
do
	mkdir bam/$a
	mkdir bam/$a/refine
	mkdir bam/$a/refine/log
	mkdir bam/$a/refine/stat/
	for bam in `cat $a.noDup.rl.bam.sorted.list`
	do
	i=`echo $bam |cut -d "." -f1 `
	echo $i
	mkdir bam/$a/refine/stat/${i}
	log="bam/$a/refine/log"
	out=${dir}/$log/qc-${i}.sh
	echo '#!/bin/bash'  > $out
	echo 'uname -a' >> $out
	echo 'watch-myself & watch_pid=$!' >> $out
	echo "java -Xmx4g -jar $GATKJAR -R $REF -knownSites $DBSNP -I ${dir}/bam/$a/refine/$bam -T CountCovariates -cov ReadGroupCovariate -cov PositionCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov CycleCovariate -recalFile ${dir}/bam/$a/refine/stat/${i}/$i.qc.csv" >> $out
	cmd="java -Xmx4g -jar $GATKCOV -recalFile ${dir}/bam/$a/refine/stat/${i}/$i.qc.csv -outputDir ${dir}/bam/$a/refine/stat/${i}/  -ignoreQ 5"
	echo $cmd >> $out
	echo 'kill -USR2 $watch_pid; wait' >> $out
	echo $cmd
	qsub -l mem=5G,time=04:: -o ${dir}/$log/qc-${i}.o -e ${dir}/$log/qc-${i}.e $out 
	done
done


########### Coverage
setting='global_setting_b37.sh'
 . $setting
#sampleList="AC4"
dir=`pwd`

for a in AC5 AC6
do
	mkdir bam/$a
	mkdir bam/$a/refine
	mkdir bam/$a/refine/log
	log="bam/$a/refine/log"
	for al in `cat $a.noDup.rl.bam.sorted.list`
	do
	bam=`basename $al`
	i=`echo $bam |cut -d "." -f1 `
	echo $i
	file="${dir}/bam/$a/refine/stat/$i.coverage"
	echo $file
	num=`grep finish ${dir}/$log/dc-${i}.o | wc -l`
	if [[ -e $file && "$num" -eq 1 ]]
	then
	echo "here $file $num"
	else
	echo "pippo"
	rm ${dir}/$log/dc-${i}.*
	if [[ "$i" =~ "GL000" ]]
	then
	echo "matched"
	out=${dir}/$log/dc-${i}.sh
	echo '#!/bin/bash'  > $out
	echo 'uname -a' >> $out
	echo 'watch-myself & watch_pid=$!' >> $out
	echo "source $myroot/$setting" >> $out
	cmd="java -Xmx2g -jar $GATKJAR -R $REF -I ${dir}/bam/$a/refine/$bam  -T DepthOfCoverage -o ${dir}/bam/$a/refine/stat/$i.coverage -ct 1 -ct 5 -ct 10 -ct 15 -ct 20 -L ${i}.1"
	echo $cmd >> $out
	echo 'kill -USR2 $watch_pid; wait' >> $out
	echo $cmd
	qsub -l mem=2G,time=02:: -o ${dir}/$log/dc-${i}.o -e ${dir}/$log/dc-${i}.e $out
	else
	echo "real chroms"
	out=${dir}/$log/dc-${i}.sh
	echo '#!/bin/bash'  > $out
	echo 'uname -a' >> $out
	echo 'watch-myself & watch_pid=$!' >> $out
	echo "source $myroot/$setting" >> $out
	cmd="java -Xmx6g -jar $GATKJAR -R $REF -I ${dir}/bam/$a/refine/$bam  -T DepthOfCoverage -o ${dir}/bam/$a/refine/stat/$i.coverage -ct 1 -ct 5 -ct 10 -ct 15 -ct 20 -L ${i}"
	echo $cmd >> $out
	echo 'kill -USR2 $watch_pid; wait' >> $out
	echo $cmd
	qsub -l mem=7G,time=20:: -o ${dir}/$log/dc-${i}.o -e ${dir}/$log/dc-${i}.e $out
	fi
	fi
	done
done

###
########### stats 
dir=`pwd`
setting='global_setting_b37.sh'
sampleList="sample.list.all.all"
. $setting
for a in AC4
do
mkdir -p ${dir}/bam/$a/refine/stat/summary
table="${dir}/bam/$a/refine/stat/summary"
rm ${table}/*
for i in {1..22} X Y MT
do
awk '{if (NR > 1) print $0}' ${dir}/bam/$a/refine/stat/${i}.coverage.sample_summary >> ${table}/summary.cov.txt
awk '{if (NR == "8") print $0}' ${dir}/bam/${a}/refine/stat/${i}/${i}.noDup.rl.sorted.fxmate_details >> ${table}/summary.mate.txt
awk -v VAR=$i '{print VAR"\t"$1}'  ${dir}/bam/$a/refine/stat/${i}.${a}.noDup.rl.flagstat > ${table}/$i-number
echo "done simple $a $i"
done
paste ${table}/*-number > ${table}/summary.mapAndpaied.txt 
~/bin/my-shuf ${dir}/bam/$a/refine/stat/1.coverage 10000 > ${table}/1-coverage.samples
echo "done shufling $a $i"
done

#######################CALLING VARIANTS
##samtool callingVariants for sample
######################################################
dir=`pwd`
setting=global_setting_b37.sh
sampleList="sample.list.all.all"
. $setting
for a in `cat $sampleList`
do
	mkdir -p bam/$a/refine/calling
	mkdir -p bam/$a/refine/calling/log
	mkdir -p ${dir}/temp
	temp="${dir}/temp"
	log="bam/$a/refine/calling/log"
	for al in `cat $a.noDup.rl.bam.sorted.list`
	do
	bam=`basename $al`
	i=`echo $bam |cut -d "." -f1 `
	echo $i
	file="${dir}/bam/$a/refine/calling/${i}.var.flt.6.vcf"
	echo $file
	num=0 #`grep finish ${dir}/$log/calSNV-${i}.o | wc -l`
	if [[ -e $file && "$num" -eq 1 ]]
	then
	echo "here $file $num"
	else
	echo "$file $num"
	if [[ "$i" =~ "GL000" ]]
	then
	echo "matched"
	#rm ${dir}/$log/calSNV-${i}.*
	out=${dir}/$log/calSNV-${i}.sh
	echo '#!/bin/bash'  > $out
	echo 'uname -a' >> $out
	echo 'watch-myself & watch_pid=$!' >> $out
	echo "source $myroot/$setting" >> $out
	#cmd="samtools mpileup -Eug -d 400 -q 1 -C 50 -6 -f $REF ${dir}/bam/$a/refine/$bam | bcftools view -bvcg - > ${dir}/bam/$a/refine/calling/${i}.var.raw.6.bcf"
	#echo $cmd >> $out
	#cmd="bcftools view ${dir}/bam/$a/refine/calling/${i}.var.raw.6.bcf | vcfutils.pl varFilter -D 400 >  ${dir}/bam/$a/refine/calling/${i}.var.flt.6.vcf"
	#echo $cmd >> $out
	#echo "java -Xmx3g -Djava.io.tmpdir=$temp -jar $GATKJAR -R $REF -T VariantAnnotator --variant ${dir}/bam/$a/refine/calling/${i}.var.flt.6.vcf  -o  ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.vcf  -all -XA MVLikelihoodRatio -XA RodRequiringAnnotation    -XA TransmissionDisequilibriumTest -XA ChromosomeCounts -XA HardyWeinberg -XA SnpEff -XA NBaseCount   -XA BaseCounts  -XA TechnologyComposition  -XA SampleList -L $i.1" >> $out
	#echo "convert2annovar.pl ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.vcf -format vcf4 -includeinfo -allallele >  ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.annovar" >> $out
	echo "summarize_annovar1.pl  ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.annovar  /ifs/scratch/c2b2/ngs_lab/sz2317/softwares/annovar/humandb/  -outfile ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.annovar.summary -ver1000g 1000g2012feb -verdbsnp 132 -genetype=refgene --buildver hg19" >> $out
	echo "perl /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/convert_annovar_vcf-all-samples.pl ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.annovar.summary.genome_summary.csv  ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.vcf" >> $out
	echo "perl /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/convert_annovar_vcf-all-samples.pl ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.annovar.summary.exome_summary.csv ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.vcf" >> $out
	echo "cp ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.annovar.summary.genome_summary.csv.vcf ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.annovar.summary.genome_summary.csv.bk.vcf">> $out
	echo "bgzip ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.annovar.summary.genome_summary.csv.vcf">> $out
	echo "tabix ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.annovar.summary.genome_summary.csv.vcf" >> $out
	echo 'kill -USR2 $watch_pid; wait' >> $out
	qsub -l mem=4G,time=2:: -o ${dir}/$log/calSNV-${i}.o -e ${dir}/$log/calSNV-${i}.e $out
	else
	#rm ${dir}/$log/calSNV-${i}.*
	out=${dir}/$log/calSNV-${i}.sh
	echo '#!/bin/bash'  > $out
	echo 'uname -a' >> $out
	echo 'watch-myself & watch_pid=$!' >> $out
	echo "source $myroot/$setting" >> $out
	#cmd="samtools mpileup -Eug -d 400 -q 1 -C 50 -6 -f $REF ${dir}/bam/$a/refine/$bam | bcftools view -bvcg - > ${dir}/bam/$a/refine/calling/${i}.var.raw.6.bcf"
	#echo $cmd >> $out
	#cmd="bcftools view ${dir}/bam/$a/refine/calling/${i}.var.raw.6.bcf | vcfutils.pl varFilter -D 400 >  ${dir}/bam/$a/refine/calling/${i}.var.flt.6.vcf"
	#echo $cmd >> $out
	#echo "java -Xmx6g -Djava.io.tmpdir=$temp -jar $GATKJAR -R $REF -T VariantAnnotator --variant ${dir}/bam/$a/refine/calling/${i}.var.flt.6.vcf  -o  ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.vcf  -all -XA MVLikelihoodRatio -XA RodRequiringAnnotation    -XA TransmissionDisequilibriumTest -XA ChromosomeCounts -XA HardyWeinberg -XA SnpEff -XA NBaseCount  -XA BaseCounts  -XA TechnologyComposition  -XA SampleList -L $i" >> $out
	#echo "convert2annovar.pl ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.vcf -format vcf4 -includeinfo -allallele >  ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.annovar" >> $out
	echo "summarize_annovar1.pl  ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.annovar  /ifs/scratch/c2b2/ngs_lab/sz2317/softwares/annovar/humandb/  -outfile ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.annovar.summary -ver1000g 1000g2012feb -verdbsnp 132 -genetype=refgene --buildver hg19" >> $out
	echo "perl /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/convert_annovar_vcf-all-samples.pl ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.annovar.summary.genome_summary.csv ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.vcf" >> $out
	echo "perl /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/convert_annovar_vcf-all-samples.pl ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.annovar.summary.exome_summary.csv ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.vcf" >> $out
	echo "cp ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.annovar.summary.genome_summary.csv.vcf ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.annovar.summary.genome_summary.csv.bk.vcf">> $out
	echo "bgzip ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.annovar.summary.genome_summary.csv.vcf">> $out
	echo "tabix ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.annovar.summary.genome_summary.csv.vcf.gz" >> $out
	echo "cp ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.annovar.summary.exome_summary.csv.vcf ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.annovar.summary.exome_summary.csv.bk.vcf">> $out
	echo "bgzip ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.annovar.summary.exome_summary.csv.vcf">> $out
	echo "tabix ${dir}/bam/$a/refine/calling/${i}.var.flt.6.ann.annovar.summary.exome_summary.csv.vcf.gz" >> $out

	echo 'kill -USR2 $watch_pid; wait' >> $out
	qsub -l mem=4G,time=5:: -o ${dir}/$log/calSNV-${i}.o -e ${dir}/$log/calSNV-${i}.e $out
	fi
	fi
	done
done

## merge vcf files) per samples 
dir=`pwd`
setting='global_setting_b37.sh'
sampleList="sample.list.all.all"
source ${dir}/$setting
for a in `cat $sampleList`
do
	mkdir -p bam/$a/refine/calling/final
	rm ${dir}/bam/$a/refine/calling/final/log.sh
	out=${dir}/bam/$a/refine/calling/final/log.sh
	#touch ${dir}/am/$a/refine/calling/final/log.sh
	#rm bam/$a/refine/calling/GL0*.var.flt.6.ann.annovar*_summary.csv.vcf.gz
	#/bin/ls bam/$a/refine/calling/*.var.flt.6.ann.annovar.summary.genome_summary.csv.vcf.gz > $dir/$a.var.flt.6.genome.vcf.gz.list
	#/bin/ls bam/$a/refine/calling/*.var.flt.6.ann.annovar.summary.exome_summary.csv.vcf.gz > $dir/$a.var.flt.6.exome.vcf.gz.list
	#awk -v var="${dir}/bam/$a/refine/calling/" '{print var$i}' $dir/var.flt.vcf.6.list.standard > $dir/$a.var.flt.6.vcf.gz.list
	#echo "vcf-concat -f $dir/$a.var.flt.6.genome.vcf.gz.list | bgzip -c > $dir/bam/$a/refine/calling/final/$a.flt.6.vcf.genome.ann.vcf.gz" >> $out
	#echo "gunzip $dir/bam/$a/refine/calling/final/$a.flt.6.vcf.genome.ann.vcf.gz" >> $out
	echo "ruby /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/filter_vcf_coding.rb -v $dir/bam/$a/refine/calling/final/$a.flt.6.vcf.genome.ann.vcf --dp4 2  --pv41 0.000001 --pv42 0.000001 --pv43 0.00001 --pv44 0.000001 " >> $out
	#echo "/ifs/data/c2b2/ngs_lab/ngs/code/NGS/utils/vcf_seperate-SNV-indel.rb $dir/bam/$a/refine/calling/final/$a.flt.6.vcf.genome.ann.filtered" >> $out
	#echo "vcf-concat -f $dir/$a.var.flt.6.exome.vcf.gz.list | bgzip -c > $dir/bam/$a/refine/calling/final/$a.flt.6.vcf.exome.ann.vcf.gz" >> $out
	#echo "gunzip $dir/bam/$a/refine/calling/final/$a.flt.6.vcf.exome.ann.vcf.gz" >> $out
	echo "ruby /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/filter_vcf_coding.rb -v $dir/bam/$a/refine/calling/final/$a.flt.6.vcf.exome.ann.vcf --dp4 2 --pv41 0.000001 --pv42 0.000001 --pv43 0.00001 --pv44 0.000001" >> $out
	#echo "ruby /ifs/data/c2b2/ngs_lab/ngs/code/NGS/utils/vcf_seperate-SNV-indel.rb $dir/bam/$a/refine/calling/final/$a.flt.6.vcf.exome.ann..filtered" >> $out
	echo "grep 1KG $dir/bam/$a/refine/calling/final/$a.flt.6.vcf.exome.ann.vcf.filtered -v| grep dbSNP -v | grep ESP5400 -v | grep functionalClass=synonymous -v > $dir/bam/$a/refine/calling/final/$a.flt.6.exome.ann.filtered.vcf" >> $out
	echo "perl /ifs/scratch/c2b2/ac_lab/pn2204/net/extractGenes.pl $dir/bam/$a/refine/calling/final/$a.flt.6.exome.ann.filtered.vcf $dir/bam/$a/refine/calling/final/$a.genes" >> $out
	chmod -X $out
	. $out
done

dir=`pwd`
setting='global_setting_b37.sh'
sampleList="sample.list.all.all"
source ${dir}/$setting
for a in AC1 Ac2 Ac3 AC4 AC5 AC6
do
$tt $dir/bam/$a/refine/calling/final/$a.flt.6.vcf.genome.ann.vcf >>stat-final2;
grep 1KG $dir/bam/$a/refine/calling/final/AC1.flt.6.vcf.genome.ann.vcf.filtered.vcf -v| grep dbSNP -v | grep ESP5400 -v | grep functionalClass=synonymous -v > $dir/bam/$a/refine/calling/final/AC1.flt.6.vcf.genome.ann.vcf.filtered.novel.vcf;
perl /ifs/scratch/c2b2/ac_lab/pn2204/net/extractGenes.pl $dir/bam/$a/refine/calling/final/AC1.flt.6.vcf.genome.ann.vcf.filtered.novel.vcf $dir/bam/$a/refine/calling/final/$a.genes;
grep 1KG ${dir}/bam/callingJoin/$a/final/$a.flt.6.vcf.gneome.ann.sorted.vcf.somaticfiltered.vcf -v| grep dbSNP -v | grep ESP5400 -v | grep functionalClass=synonymous -v > ${dir}/bam/callingJoin/${a}/final/${a}.flt.6.vcf.genome.ann.sorted.vcf.somaticfiltered.Known.vcf" >> $out
	echo "/ifs/scratch/c2b2/ac_lab/pn2204/net/script/extractGenesSomatic.pl ${dir}/bam/callingJoin/$a/final/${a}.flt.6.vcf.genome.ann.sorted.vcf.somaticfiltered.Known.vcf $dir/bam/callingJoin/$a/final/$a.genes" >> $out


done
rm $dir/bam/$a/refine/calling/*.var.flt.6.ann.annovar
rm $dir/bam/$a/refine/calling/*.var.flt.6.ann.annovar.summary.variant_function
rm $dir/bam/$a/refine/calling/*.var.flt.6.ann.annovar.summary.invalid*
rm $dir/bam/$a/refine/calling/*.var.flt.6.ann.annovar.summary.hg19*
rm $dir/bam/$a/refine/calling/*.var.flt.6.vcf
rm $dir/bam/$a/refine/calling/*var.raw.6.bcf
rm *.var.flt.6.ann.annovar
rm *.var.flt.6.ann.annovar.summary.variant_function
rm *.var.flt.6.ann.annovar.summary.invalid*
rm *.var.flt.6.ann.annovar.summary.hg19*
rm *.var.flt.6.vcf
rm *var.raw.6.bcf
rm *var.flt.6.vcf.idx
rm *var.flt.6.ann.annovar.summary.exonic_variant_function
rm *var.flt.6.ann.annovar.summary.log
rm *joint.var.flt.6.ann.vcf*






####################################
##concordance
##########################
 java -Xmx6g -Djava.io.tmpdir=$TEMP -jar $GATKJAR -R $REF -T VariantEval -eval ${dir}/bam/AC5/refine/calling/22.flt.6.vcf  --evalModule --known_names ${dir}/bam/AC2/refine/calling/22.flt.6.ann.vcf -o 22.eval.gatkreport 

java -Xmx6g  -jar $GATKJAR -R $REF -T VariantEval -eval:AC5 ${dir}/bam/AC5/refine/calling/22.var.flt.6.vcf --eval:AC6 ${dir}/bam/AC6/refine/calling/22.var.flt.6.vcf --eval --comp:AC1 ${dir}/bam/AC1/refine/calling/22.var.flt.6.ann.vcf --comp:AC2 ${dir}/bam/AC2/refine/calling/22.var.flt.6.ann.vcf  --comp:AC3 ${dir}/bam/AC3/refine/calling/22.var.flt.6.ann.vcf --comp:AC4 ${dir}/bam/AC4/refine/calling/22.var.flt.6.ann.vcf -o 22.eval.gatkreport -L 22 -D $DBSNP
################################
##samtool callingVariants for paiered sample
###############################################
dir=`pwd`
setting=global_setting_b37.sh
sampleList="sample.list.all.all"
chrList="chrom.list"
. $setting
for a in AC2
do
	for b in AC3
	do
	if [ "$a" = "$b" ];
	then
	echo "same"
	else
	mkdir -p bam/callingJoin/$a-$b
	mkdir -p bam/callingJoin/$a-$b/log
	log="bam/callingJoin/$a-$b/log"
	for al in `cat AC1.noDup.rl.bam.sorted.list`
	do
	bam=`basename $al`
	i=`echo $bam |cut -d "." -f1 `
	echo $i
	file="${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.vcf"
	echo $file
	num=`grep finish ${dir}/$log/calSNV-${a}-${b}-${i}.o | wc -l`
	if [[ -e $file && "$num" -eq 1 ]]
	then
	echo "here $file $num"
	else
	echo "$file $num"
	if [[ "$i" =~ "GL000" ]]
	then
	echo "matched"
	#rm ${dir}/$log/calSNV-${a}-${b}-${i}.*
	#out=${dir}/$log/calSNV-${a}-${b}-${i}.2.sh
	#echo '#!/bin/bash'  > $out
	#echo 'uname -a' >> $out
	#echo 'watch-myself & watch_pid=$!' >> $out
	#echo "source $myroot/$setting" >> $out
	#cmd="samtools mpileup -DSu -d 400 -6 -q 1 -C 50 -f $REF ${dir}/bam/$a/refine/$bam ${dir}/bam/$b/refine/$bam | bcftools view -bvcgT pair - >  ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.raw.6.bcf"
	#echo $cmd >> $out
	#cmd="bcftools view ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.raw.6.bcf | vcfutils.pl varFilter -D 400 >  ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.vcf"
	#echo $cmd >> $out
	#echo "java -Xmx3g -Djava.io.tmpdir=$TMPDIR -jar $GATKJAR -R $REF -T VariantAnnotator --variant ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.vcf  -o  ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.vcf -all -XA MVLikelihoodRatio -XA RodRequiringAnnotation    -XA TransmissionDisequilibriumTest -XA ChromosomeCounts -XA HardyWeinberg -XA SnpEff -XA NBaseCount  -XA BaseCounts  -XA TechnologyComposition  -XA SampleList -L $i.1" >> $out
	#echo "convert2annovar.pl ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.vcf -format vcf4 -includeinfo -allallele > ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar" >> $out
	#echo "summarize_annovar1.pl  ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar  /ifs/scratch/c2b2/ngs_lab/sz2317/softwares/annovar/humandb/  -outfile ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary -ver1000g 1000g2012feb -verdbsnp 132 -genetype=refgene --buildver hg19" >> $out
	#echo "perl /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/convert_annovar_vcf-all-samples.pl ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary.genome_summary.csv ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.vcf" >> $out
	#echo "perl /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/convert_annovar_vcf-all-samples.pl ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary.exome_summary.csv ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.vcf" >> $out
	#echo "cp ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary.genome_summary.csv.vcf ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary.genome_summary.bk.csv.vcf">> $out
	#echo "bgzip -f ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary.genome_summary.csv.vcf">> $out
	#echo "tabix -f ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary.genome_summary.csv.vcf.gz" >> $out
	#echo "cp ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary.exome_summary.csv.vcf ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary.exome_summary.bk.csv.vcf">> $out
	#echo "bgzip -f ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary.exome_summary.csv.vcf">> $out
	#echo "tabix -f ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary.exome_summary.csv.vcf.gz" >> $out
	#echo 'kill -USR2 $watch_pid; wait' >> $out
	#qsub -l mem=4G,time=3:: -o ${dir}/$log/calSNV-${a}-${b}-${i}.o -e ${dir}/$log/calSNV-${a}-${b}-${i}.e $out
	else
	rm ${dir}/$log/calSNV-${a}-${b}-${i}.*
	out=${dir}/$log/calSNV-${a}-${b}-${i}.sh
	echo '#!/bin/bash'  > $out
	echo 'uname -a' >> $out
	echo 'watch-myself & watch_pid=$!' >> $out
	echo "source $dir/$setting" >> $out
	echo "samtools mpileup -DSu -d 400 -6 -q 1 -C 50 -f $REF ${dir}/bam/$a/refine/$bam ${dir}/bam/$b/refine/$bam | bcftools view -bvcgT pair - >  ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.raw.6.bcf" >> $out
	echo "bcftools view ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.raw.6.bcf | vcfutils.pl varFilter -D 400 >  ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.vcf" >> $out
	echo "java -Xmx9g -Djava.io.tmpdir=$TMPDIR -jar $GATKJAR -R $REF -T VariantAnnotator --variant ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.vcf  -o  ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.vcf -all -XA MVLikelihoodRatio -XA RodRequiringAnnotation    -XA TransmissionDisequilibriumTest -XA ChromosomeCounts -XA HardyWeinberg -XA SnpEff -XA NBaseCount  -XA BaseCounts  -XA TechnologyComposition  -XA SampleList -L ${i} " >> $out
	echo "convert2annovar.pl ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.vcf -format vcf4 -includeinfo -allallele > ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar " >> $out
	echo "summarize_annovar1.pl  ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar  /ifs/scratch/c2b2/ngs_lab/sz2317/softwares/annovar/humandb/  -outfile ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary -ver1000g 1000g2012feb -verdbsnp 132 -genetype=refgene --buildver hg19" >> $out
	echo "rm $dir/bam/callingJoin/$a-$b/*.var.flt.6.ann.annovar $dir/bam/callingJoin/$a-$b//*.var.flt.6.ann.annovar.summary.variant_function $dir/bam/callingJoin/$a-$b/*.var.flt.6.ann.annovar.summary.invalid* $dir/bam/callingJoin/$a-$b/*.var.flt.6.ann.annovar.summary.hg19* $dir/bam/callingJoin/$a-$b/*.var.flt.6.vcf $dir/bam/callingJoin/$a-$b/*var.raw.6.bcf " >>$out
	echo "perl /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/convert_annovar_vcf-all-samples.pl ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary.genome_summary.csv ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.vcf" >> $out
	echo "perl /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/convert_annovar_vcf-all-samples.pl ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary.exome_summary.csv ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.vcf" >> $out
	#echo "/ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/annotate_joint_vcf_dp4.rb --vcf  ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary.genome_summary.csv.vcf --s1 $i.pileup --s2 $2.pileup "
	#echo "cp ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary.genome_summary.csv.vcf ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary.genome_summary.bk.csv.vcf">> $out
	echo "bgzip -f ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary.genome_summary.csv.vcf">> $out
	echo "tabix -f ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary.genome_summary.csv.vcf.gz" >> $out
	#echo "cp ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary.exome_summary.csv.vcf ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary.exome_summary.bk.csv.vcf">> $out
	echo "bgzip -f ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary.exome_summary.csv.vcf">> $out
	echo "tabix -f ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary.exome_summary.csv.vcf.gz" >> $out
	echo 'kill -USR2 $watch_pid; wait' >> $out
	qsub -l mem=10G,time=15:: -o ${dir}/$log/calSNV-${a}-${b}-${i}.o -e ${dir}/$log/calSNV-${a}-${b}-${i}.e $out
	fi
	fi
	done
	fi
	done
done

#####################
dir=`pwd`
for a in AC1 AC2 AC3 AC4 AC5 AC6
do
ruby /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/transition-transversion_vcf.rb  $dir/bam/$a/refine/calling/final/${a}.flt.6.vcf.genome.ann.vcf > $dir/bam/$a/refine/calling/final/stat-final-2
done

## merge vcf files per pair
###################################################################

dir=`pwd`
setting='global_setting_b37.sh'
sampleList="join.call.sample"
source ${dir}/$setting
for al in `cat $sampleList`
do
	a=`basename $al`
	mkdir -p ${dir}/bam/callingJoin/$a/final
	rm ${dir}/bam/callingJoin/$a/final/log.sh
	out=${dir}/bam/callingJoin/$a/final/log.2.sh
	#touch ${dir}/bam/callingJoin/$a/final/log.sh
	#rm ${dir}/bam/callingJoin/$a/GL0*.var.flt.6.ann.annovar*_summary.csv.vcf.gz
	#/bin/ls ${dir}/bam/callingJoin/$a/*.var.flt.6.ann.annovar.summary.genome_summary.csv.vcf.gz > $dir/$a.var.flt.6.genome.vcf.gz.list
	#/bin/ls ${dir}/bam/callingJoin/$a/*var.flt.6.ann.annovar.summary.exome_summary.csv.vcf.gz > $dir/$a.var.flt.6.exome.vcf.gz.list
	#awk -v var="${dir}/bam/$a/refine/calling/" '{print var$i}' $dir/var.flt.vcf.6.list.standard > $dir/$a.var.flt.6.vcf.gz.list
	#echo "vcf-concat -f $dir/$a.var.flt.6.genome.vcf.gz.list | bgzip -c > ${dir}/bam/callingJoin/$a/final/$a.flt.6.vcf.genome.ann.vcf.gz" >> $out
	#echo "vcf-sort ${dir}/bam/callingJoin/$a/final/$a.flt.6.vcf.genome.ann.vcf.gz > ${dir}/bam/callingJoin/$a/final/$a.flt.6.vcf.genome.ann.sorted.vcf" >>$out
	#echo "ruby /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/filter_vcf_somatic.rb  -v ${dir}/bam/callingJoin/$a/final/$a.flt.6.vcf.genome.ann.sorted.vcf --dp4 2 --pv41 0.000001 --pv42 0.000001 --pv43 0.00001 --pv44 0.000001 --normal 1" >> $out
	#echo "ruby /ifs/data/c2b2/ngs_lab/ngs/code/NGS/utils/vcf_seperate-SNV-indel.rb $dir/bam/callingJoin/$a/final/$a.flt.6.vcf.genome.ann.sorted.vcf.somaticfiltered.vcf" >> $out
	echo "grep 1KG ${dir}/bam/callingJoin/$a/final/$a.flt.6.vcf.genome.ann.sorted.vcf.somaticfiltered.vcf -v| grep dbSNP -v | grep ESP5400 -v | grep functionalClass=synonymous -v > ${dir}/bam/callingJoin/${a}/final/${a}.flt.6.vcf.genome.ann.sorted.vcf.somaticfiltered.Known.vcf" >> $out
	echo "/ifs/scratch/c2b2/ac_lab/pn2204/net/script/extractGenesSomatic.pl ${dir}/bam/callingJoin/$a/final/${a}.flt.6.vcf.genome.ann.sorted.vcf.somaticfiltered.Known.vcf $dir/bam/callingJoin/$a/final/$a.genes.allgenome" >> $out

	##exome
	#echo "vcf-concat -f $dir/$a.var.flt.6.exome.vcf.gz.list | bgzip -c > $dir/bam/callingJoin/$a/final/$a.flt.6.vcf.exome.ann.vcf.gz" >> $out
	#echo "vcf-sort ${dir}/bam/callingJoin/$a/final/$a.flt.6.vcf.exome.ann.vcf.gz > ${dir}/bam/callingJoin/$a/final/$a.flt.6.vcf.exome.ann.sorted.vcf" >> $out
	#echo "ruby /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/filter_vcf_somatic.rb  -v ${dir}/bam/callingJoin/$a/final/$a.flt.6.vcf.exome.ann.sorted.vcf --dp4 2 --pv41 0.000001 --pv42 0.000001 --pv43 0.00001 --pv44 0.000001 --normal 1" >> $out
	#echo "grep 1KG ${dir}/bam/callingJoin/$a/final/$a.flt.6.vcf.exome.ann.sorted.vcf.somaticfiltered.vcf -v| grep dbSNP -v | grep ESP5400 -v | grep functionalClass=synonymous -v > ${dir}/bam/callingJoin/${a}/final/${a}.flt.6.vcf.exome.ann.sorted.vcf.somaticfiltered.Known.vcf" >> $out
	#echo "/ifs/scratch/c2b2/ac_lab/pn2204/net/script/extractGenesSomatic.pl ${dir}/bam/callingJoin/$a/final/${a}.flt.6.vcf.exome.ann.sorted.vcf.somaticfiltered.Known.vcf $dir/bam/callingJoin/$a/final/$a.genes" >> $out
	echo "ruby /ifs/data/c2b2/ngs_lab/ngs/code/NGS/utils/vcf_seperate-SNV-indel.rb $dir/bam/callingJoin/$a/final/$a.flt.6.vcf.exome.ann.sorted.vcf.somaticfiltered.vcf" >> $out
	chmod -X $out
	. $out
done

### 

--pv41 0.000001 --pv42 0.00000001 --pv43 0.0001 --pv44 0.0001 --clr 40

#################
## FILTERING PARAMITERS
##############################
#Include only sites with all Non-Reference Allele Counts within the specified range > 2

vcf-stats AC1-AC2-jc.vcf.gz
vcftools --gzvcf AC1-AC2-jc.vcf.gz --site-quality --out quality
vcf-
vcftools --gzvcf AC1-AC2-jc.vcf.gz  --minQ 20 --non-ref-ac 3 --mac 0 --min-meanDP 3 --out AC1-AC2-jc.vcf.fiter.3  --recode-INFO CLR --recode-INFO CLR --recode

#######filtering somatic
##filter 1)  stradb=0.0001 ,  mapQ=0.00000001,  baseQ=0.000001
ruby /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/filter_vcf_somatic.rb --vcf call.joint.var.flt.6.vcf  --dp4 2 --pv41 0.000001 --pv42 0.00000001 --pv43 0.0001 --pv44 0.0000001 --clr 20 (old)
##more stringent filter
ruby /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/filter_vcf_somatic.rb --vcf call.joint.var.flt.6.vcf  --dp4 2 --pv41 0.0001 --pv42 0.0001 --pv43 0.0001 --pv44 0.0001 --clr 20
##
ruby /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/filter_vcf_somatic.rb --vcf call.joint.var.flt.6.vcf  --dp4 0 --pv41 0.000001 --pv42 0.00000001 --pv43 0.0001 --pv44 0.0001 --clr 40

ruby /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/filter_annovar_coding.rb --annovar FILE.txt --dp4 2 --pv41 0.0001 --pv42 0.0001 --pv43 0.0001 --pv44 0.0001


SNVs:
Consensus score >=20
Quality score >=20
A minimum of 3 reads supporting the variant

Indels:
Consensus score >=20
Quality score >=50 (FQ)
A minimum of 3 reads supporting the variant


####################################################
#######calling with VARSCAN
################################################
#from bam file procude mpilep file with samtools
dir=`pwd`
setting='global_setting_b37.sh'
sampleList="sample.list.all"
source ./$setting
mkdir -p snv/varscan/
mkdir -p snv/varscan/log
logs='snv/varscan/log'

for a in `cat $sampleList`
do
	for al in `cat $a.noDup.rl.bam.sorted.list`
	do
	bam=`basename $al`
	echo $bam
	i=`echo $bam |cut -d "." -f1 `
	echo $i
	out=${dir}/$logs/mpileup-$a-${i}.sh
	echo '#!/bin/bash'  >> $out
	echo 'uname -a' >> $out
	echo "source $myroot/$setting" >> $out
	echo 'watch-myself & watch_pid=$!' >> $out
	cmd="samtools mpileup -q 1 -f $REF ${dir}/$al > ${dir}/snv/varscan/$i.$a.pileup"
	echo $cmd >> $out
	cmd='kill -USR2 $watch_pid; wait' 
	echo $cmd >> $out
	qsub -l mem=10G,time=50:: -o ${dir}/$logs/mpilep${i}.o -e ${dir}/$logs/mpilep${i}.e $out
	done
done
##calling variants
# we called alll the variance with 
/bin/ls /snv/varscan/*.pileup > bam.pileup.list
dir=`pwd`
caselist='case.list'
controlist='control.list'
chromlist='chrom.list'
setting='global_setting_b37.sh'
 . $setting
mkdir -p snv/varscan/calling
mkdir -p snv/varscan/calling/log
logs='snv/varscan/calling/log'
for a in `cat $caselist`; 
do
	for b in `cat $controlist`; 
	do
	for i in `cat $chromlist`
	do
	out=${dir}/$logs/calling-${i}.sh
	echo '#!/bin/bash'  >> $out
	echo 'uname -a' >> $out
	echo "source $myroot/$setting" >> $out
	echo 'watch-myself & watch_pid=$!' >> $out
	cmd="java -jar $varscan somatic ${dir}/snv/varscan/$i.$b.pileup ${dir}/snv/varscan/$i.$a.pileup ${dir}/snv/varscan/calling/$i.basename --strand-filter 1  --min-coverage 4 --min-coverage-normal 2 --min-coverage-tumor 2 --min_var_freq 0.20 --p-value 1.0e-01 --somatic-p-value 1.0e-04 --tumor-purity 0.95"
	echo $cmd >> $out
	echo $cmd
	cmd="java -jar $varscan processSomatic ${dir}/snv/varscan/calling/$i.basename.snp " 
	echo $cmd >> $out
	echo $cmd
	cmd='kill -USR2 $watch_pid; wait' 
	echo $cmd >> $out
	qsub -l mem=10G,time=50:: -o ${dir}/$logs/calling-${i}.o -e ${dir}/$logs/calling-${i}.e $out
done
done
done
# we called alll the variance with stringent cut off
dir=`pwd`
caselist='case.list'
controlist='control.list'
chromlist='chrom.list'
setting='global_setting_b37.sh'
 . $setting
mkdir -p snv/varscan/calling
mkdir -p snv/varscan/calling/log
logs='snv/varscan/calling/log'
for a in `cat $caselist`; 
do
	for b in `cat $controlist`; 
	do
	for i in `cat $chromlist`
	do
	out=${dir}/$logs/calling2-${i}.sh
	echo '#!/bin/bash'  >> $out
	echo 'uname -a' >> $out
	echo "source $myroot/$setting" >> $out
	echo 'watch-myself & watch_pid=$!' >> $out
	cmd="java -jar $varscan somatic ${dir}/snv/varscan/$i.$b.pileup ${dir}/snv/varscan/$i.$a.pileup ${dir}/snv/varscan/calling/$i.basename.str --strand-filter 1  --min-coverage 10 --min-coverage-normal 8 --min-coverage-tumor 8 --min_var_freq 0.20 --p-value 1.0e-01 --somatic-p-value 1.0e-04 --tumor-purity 0.95"
	echo $cmd >> $out
	echo $cmd
	cmd="java -jar $varscan processSomatic ${dir}/snv/varscan/calling/$i.basename.str.snp " 
	echo $cmd >> $out
	echo $cmd
	cmd='kill -USR2 $watch_pid; wait' 
	echo $cmd >> $out
	qsub -l mem=10G,time=50:: -o ${dir}/$logs/calling2-${i}.o -e ${dir}/$logs/calling2-${i}.e $out
done
done
done

##Annotated
dir=`pwd`
caselist='case.list'
controlist='control.list'
chromlist='chrom.list'
setting='global_setting_b37.sh'
 . $setting
mkdir -p ${dir}/snv/varscan/calling/final
mkdir -p ${dir}/snv/varscan/calling/final/log
log="snv/varscan/calling/final/log"
for i in `cat $chromlist`
do
awk '{if (NR>1) print $0 }' ${dir}/snv/varscan/calling/$i.basename.snp.Somatic.hc >> ${dir}/snv/varscan/calling/final/all.snp.Somatic.hc
awk '{if (NR>1) print $0 }' ${dir}/snv/varscan/calling/$i.basename.str.snp.Somatic.hc >> ${dir}/snv/varscan/calling/final/all.str.snp.Somatic.hc
awk '{if (NR>1) print $0 }' ${dir}/snv/varscan/calling/$i.basename.indel >> ${dir}/snv/varscan/calling/final/all.indel
awk '{if (NR>1) print $0 }' ${dir}/snv/varscan/calling/$i.basename.str.indel >> ${dir}/snv/varscan/calling/final/all.str.indel
awk '{if (NR>1) print $0 }' ${dir}/snv/varscan/calling/$i.basename.snp.Germline >> ${dir}/snv/varscan/calling/final/all.snp.Germline
awk '{if (NR>1) print $0 }' ${dir}/snv/varscan/calling/$i.basename.snp.LOH >> ${dir}/snv/varscan/calling/final/all.basename.snp.LOH
awk '{if (NR>1) print $0 }' ${dir}/snv/varscan/calling/$i.basename.str.snp.Germline >> ${dir}/snv/varscan/calling/final/all.str.snp.Germline
awk '{if (NR>1) print $0 }' ${dir}/snv/varscan/calling/$i.basename.str.snp.LOH >> ${dir}/snv/varscan/calling/final/all.str.snp.LOH
done 

dir=`pwd`
caselist='case.list'
controlist='control.list'
chromlist='chrom.list'
setting='global_setting_b37.sh'
 . $setting
mkdir -p ${dir}/snv/varscan/calling/final
mkdir -p ${dir}/snv/varscan/calling/final/log
log="snv/varscan/calling/final/log"
/bin/ls snv/varscan/calling/final/all* > all.varscan.list
for al in `cat all.varscan.list`
do
i=`basename $al`
awk '{print $1"\t"$2"\t"$2"\t"$3"\t"$4}' ${dir}/snv/varscan/calling/final/$i > $i.annovar
done

dir=`pwd`
caselist='case.list'
controlist='control.list'
chromlist='chrom.list'
setting='global_setting_b37.sh'
 . $setting

#/bin/ls all*.annovar> list
for i in `cat list`
do
out=${dir}/$log/an-$i.sh
echo '#!/bin/bash' >> $out
echo 'uname -a' >> $out
echo "source $myroot/$setting" >> $out
echo 'watch-myself & watch_pid=$!' >> $out
echo "summarize_annovar.pl --outfile ${dir}/snv/varscan/calling/final/$i.annovar.summary ${dir}/snv/varscan/calling/final/$i /ifs/scratch/c2b2/ngs_lab/sz2317/softwares/annovar/humandb/  -ver1000g 1000g2012feb -verdbsnp 132 -genetype=refgene --buildver hg19" >> $out
echo 'kill -USR2 $watch_pid; wait' >> $out
qsub -l mem=2G,time=5:: -o ${dir}/$log/an-$i.o -e ${dir}/$log/an-$i.e $out
done

######################################
#########Exon data
##############################################################################
dir=`pwd`
setting='global_setting_b37.sh'

 . $setting

for a in AC7 AC8 AC9
do
	mkdir -p bam/$a
	rm bam/$a/*
	out=${dir}/bam/$a/log.sh
	touch ${dir}/bam/$a/log.sh
	cp /ifs/scratch/c2b2/ngs_lab/ngs/Projects/Exome-seq/120309_ANDREA_GABRIELLE_3_HUMAN_EXOME_60X_PE_HISEQ/120419_SN828_0129_AD0R2DACXX/release/$a.indel.vcf  bam/$a
	cp /ifs/scratch/c2b2/ngs_lab/ngs/Projects/Exome-seq/120309_ANDREA_GABRIELLE_3_HUMAN_EXOME_60X_PE_HISEQ/120419_SN828_0129_AD0R2DACXX/release/$a.SNV.vcf bam/$a
	echo "bgzip bam/$a/$a.SNV.vcf">> $out
	echo "tabix bam/$a/$a.SNV.vcf.gz" >> $out
	echo "bgzip bam/$a/$a.indel.vcf">> $out
	echo "tabix bam/$a/$a.indel.vcf.gz" >> $out
	echo "done copy $a"
	echo "vcf-concat bam/$a/$a.SNV.vcf.gz bam/$a/$a.indel.vcf.gz| gzip -c > bam/$a/temp.gz " >> $out
	echo "vcf-sort bam/$a/temp.gz >  bam/$a/$a.vcf" >> $out
	echo "done concate $a"
	echo "rm bam/$a/temp.gz " >> $out
	#echo "gunzip bam/$a/$a.vcf.gz" >> $out
	echo "grep "1KG" bam/$a/$a.vcf -v| grep dbSNP -v | grep EVS -v| grep functionalClass=silent -v |grep intron -v | grep refseq.functionalClass=noncoding -v |grep PASS| grep 0/1 > bam/$a/$a.filtered.woLOH.list" >> $out
	echo "perl /ifs/scratch/c2b2/ac_lab/pn2204/net/script/extractGenesAnnotationGATK.pl bam/$a/$a.filtered.vcf bam/$a/$a.genes" >> $out
	echo "done  $a"
	chmod -X $out
	. $out
done

for a in AC7 AC8 AC9 
do 
grep "1KG" $a/$a.vcf -v | grep dbSNP -v | grep EVS -v| grep functionalClass=silent -v |grep intron -v | grep refseq.functionalClass=noncoding -v | grep utr5 -v | grep utr3 -v |grep refseq.functionalClass=inframe -v | grep PASS | grep 0/1 > $a/$a.filtered.woLOH.list; 
perl ../script/extractGenesAnnotationGATK.pl $a/$a.filtered.woLOH.list $a/$a.gene.new; 
done

java -Xmx6g  -jar $GATKJAR -R $REF -T VariantEval -eval:AC5 ${dir}/bam/AC5/refine/calling/final/AC5.flt.6.vcf.genome.ann.vcf --eval:AC6 ${dir}/bam/AC6/refine/calling/final/AC5.flt.6.vcf.genome.ann.vcf --comp:AC7 ${dir}/bam/AC1/refine/calling/22.var.flt.6.ann.vcf --comp:AC8 ${dir}/bam/AC2/refine/calling/22.var.flt.6.ann.vcf  --comp:AC9 ${dir}/bam/AC3/refine/calling/22.var.flt.6.ann.vcf --comp:AC4 ${dir}/bam/AC4/refine/calling/22.var.flt.6.ann.vcf -o 22.eval.gatkreport -L 22 -D $DBSNP

mkdir -p clSimilarity
out=${dir}/bam/$a/log.sh
touch ${dir}/bam/$a/log.sh
echo  "java -Xmx6g -jar $GATKJAR -R $REF -T VariantEval -eval:AC5 ${dir}/bam/AC5/refine/calling/final/AC5.flt.6.vcf.genome.ann.vcf --eval:AC6 ${dir}/bam/AC6/refine/calling/final/AC6.flt.6.vcf.genome.ann.vcf --comp:AC7 ${dir}/bam/AC7/AC7.vcf --comp:AC8 ${dir}/bam/AC8/AC8.vcf  --comp:AC9 ${dir}/bam/AC9/AC9.vcf  -o clSimilarity/CellLineEvalu.eval.gatkreport -D $DBSNP" >>$out


################################
##samtool callingVariants for paiered exome sample
###############################################
dir=`pwd`
setting=global_setting_b37.sh
dirG="/ifs/scratch/c2b2/ngs_lab/ngs/Projects/Exome-seq/120309_ANDREA_GABRIELLE_3_HUMAN_EXOME_60X_PE_HISEQ/120419_SN828_0129_AD0R2DACXX/release/"
sampleList="sample.list.all.all"
chrList="chrom.list"
. $setting
for a in AC9
do
	for b in AC7
	do
	if [ "$a" = "$b" ];
	then
	echo "same"
	else
	mkdir -p bam/callingJoin/$a-$b
	mkdir -p bam/callingJoin/$a-$b/log
	log="bam/callingJoin/$a-$b/log"
	file="${dir}/bam/callingJoin/$a-$b/$a-$b.joint.var.flt.6.ann.annovar.summary.exome_summary.csv.vcf.gz"
	echo $file
	num=0 #`grep finish ${dir}/$log/calSNV-${a}-${b}-${i}.o | wc -l`
	if [[ -e $file && "$num" -eq 1 ]]
	then
	echo "here $file $num"
	else
	echo "$file $num"
	rm ${dir}/$log/calSNV-${a}-${b}-${i}.*
	out=${dir}/$log/calSNV-${a}-${b}-${i}.sh
	echo '#!/bin/bash'  > $out
	echo 'uname -a' >> $out
	echo 'watch-myself & watch_pid=$!' >> $out
	echo "source $myroot/$setting" >> $out
	echo "samtools mpileup -DSu -d 400 -q 1 -C 50 -f ${Exome} ${dirG}/${a}.bam ${dirG}/${b}.bam | bcftools view -bvcgT pair - >  ${dir}/bam/callingJoin/$a-$b/$a-$b.joint.var.raw.6.bcf" >> $out
	echo "bcftools view ${dir}/bam/callingJoin/$a-$b/$a-$b.joint.var.raw.6.bcf | vcfutils.pl varFilter -D 400 >  ${dir}/bam/callingJoin/$a-$b/$a-$b.joint.var.flt.6.vcf" >> $out
	echo "java -Xmx4g -Djava.io.tmpdir=$TMPDIR -jar $GATKJAR -R $REF -T VariantAnnotator --variant ${dir}/bam/callingJoin/$a-$b/$a-$b.joint.var.flt.6.vcf  -o  ${dir}/bam/callingJoin/$a-$b/$a-$b.joint.var.flt.6.ann.vcf -all -XA MVLikelihoodRatio -XA RodRequiringAnnotation    -XA TransmissionDisequilibriumTest -XA ChromosomeCounts -XA HardyWeinberg -XA SnpEff -XA NBaseCount  -XA BaseCounts  -XA TechnologyComposition  -XA SampleList -L $Exome " >> $out
	echo "convert2annovar.pl ${dir}/bam/callingJoin/$a-$b/$a-$b.joint.var.flt.6.ann.vcf -format vcf4 -includeinfo -allallele > ${dir}/bam/callingJoin/$a-$b/$a-$b.joint.var.flt.6.ann.annovar" >> $out
	echo "summarize_annovar1.pl  ${dir}/bam/callingJoin/$a-$b/$a-$b.joint.var.flt.6.ann.annovar  /ifs/scratch/c2b2/ngs_lab/sz2317/softwares/annovar/humandb/  -outfile ${dir}/bam/callingJoin/$a-$b/$a-$b.joint.var.flt.6.ann.annovar.summary -ver1000g 1000g2012feb -verdbsnp 132 -genetype=refgene --buildver hg19" >> $out
	echo "rm $dir/bam/callingJoin/$a-$b/*.var.flt.6.ann.annovar $dir/bam/callingJoin/$a-$b/*.var.flt.6.ann.annovar.summary.variant_function $dir/bam/callingJoin/$a-$b/*.var.flt.6.ann.annovar.summary.invalid* $dir/bam/callingJoin/$a-$b/*.var.flt.6.ann.annovar.summary.hg19* $dir/bam/callingJoin/$a-$b/*.var.flt.6.vcf $dir/bam/callingJoin/$a-$b/*var.raw.6.bcf " >>$out
	echo "perl /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/convert_annovar_vcf-all-samples.pl ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.annovar.summary.exome_summary.csv ${dir}/bam/callingJoin/$a-$b/${i}.joint.var.flt.6.ann.vcf" >> $out
	echo "bgzip -f ${dir}/bam/callingJoin/$a-$b/$a-$b.joint.var.flt.6.ann.annovar.summary.exome_summary.csv.vcf">> $out
	echo "tabix -f ${dir}/bam/callingJoin/$a-$b/$a-$b.joint.var.flt.6.ann.annovar.summary.exome_summary.csv.vcf.gz" >> $out
	echo 'kill -USR2 $watch_pid; wait' >> $out
	qsub -l mem=5G,time=20:: -o ${dir}/$log/calSNV-${a}-${b}-${i}.o -e ${dir}/$log/calSNV-${a}-${b}-${i}.e $out
	fi
	fi
	done
done

