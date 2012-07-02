##calling variant from AC1 (normal) and AC2 (liver met)

##level 1 fasta q file




##level 2 bam files



##level 3 calling variant

##validate the calling by picard

REF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/human_g1k_v37.fasta"
dir=`pwd`
script="/ifs/scratch/c2b2/ac_lab/pn2204/net/script/picard_validateSam.sh"
Inputdir="/ifs/scratch/c2b2/ngs_lab/ngs/Projects/WG-seq/120309_ANDREA_GABRIELLE_6_HUMAN_GENOME_30X_PE_HISEQ/SAMPLES/AC1/BAM"
sampleList="sample.list"
bamlist="bam.name.list"
mkdir bam
for a in `cat $sampleList`
do
	mkdir bam/$a
	mkdir bam/$a/validBam
	mkdir bam/$a/validBam/log
	log="bam/$a/validBam/log"
	for bam in `cat $bamlist`
	do
	i=`echo $bam |cut -d "." -f1 `
	out=${dir}/$log/vaLPic-${i}.sh
	echo '#!/bin/bash'  > $out
	echo 'uname -a' >> $out
	cmd="$script -i ${Inputdir}/$bam -o ${dir}/bam/$a/validBam/$i.val -r $REF "
	echo $cmd >> $out
	echo $cmd
	echo $out
	qsub -l mem=5G,time=50:: -o ${dir}/$log/vaLPic${i}.o -e ${dir}/$log/vaLPic${i}.e $out
	done
done
##remove PCR duplicate

samtool="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/samtools"
Inputdir="/ifs/scratch/c2b2/ngs_lab/ngs/Projects/WG-seq/120309_ANDREA_GABRIELLE_6_HUMAN_GENOME_30X_PE_HISEQ/SAMPLES/AC1/BAM"
sampleList="sample.list.AC1"
chrmlist="chrom.list"
bamlist="bam.name.list"
mkdir bam
dir=`pwd`
#from bam file elimina PCR duplicate read with samtools for each sample for each chromosome
for a in `cat $sampleList`
do
	mkdir bam/$a
	mkdir bam/$a/log
	for bam in `cat $bamlist`
	do
	echo $bam
	i=`echo $bam $|cut -d "." -f1 `
	echo $i
	out=${dir}/bam/$a/log/removeDupl${i}.sh
	echo '#!/bin/bash'  > $out
	echo 'uname -a' >> $out
	cmd="$samtool rmdup  ${Inputdir}/$bam ${dir}/bam/$a/${i}.noDup.bam "
	echo $cmd >> $out
	echo $cmd
	qsub -l mem=4G,time=5:: -o ${dir}/bam/$a/log/removeDupl${i}.o -e ${dir}/bam/$a/log/removeDupl${i}.e $out
	done
done
##called is the samtool to reallined
ls bam/AC1/*bam | cut -d "/" -f3 > bam.noDup.list
REF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/human_g1k_v37.fasta"
samtool="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/samtools"
sampleList="sample.list.AC1"
bamlist="bam.noDup.list"

dir=`pwd`
for a in `cat $sampleList`
do
	mkdir bam/$a
	mkdir bam/$a/refine
	mkdir bam/$a/refine/log
	rm ${dir}/bam/$a/refine/log/* ##attention it is a rm option
	for bam in `cat $bamlist`
	do
	i=`echo $bam |cut -d "." -f1 `
	echo $i
	log="bam/$a/refine/log"
	out=${dir}/$log/calmd-${i}.sh
	echo '#!/bin/bash'  > $out
	echo 'uname -a' >> $out
	echo 'watch-myself &' >> $out
	cmd="$samtool calmd -rEA ${dir}/bam/${a}/${bam} $REF > ${dir}/bam/$a/refine/${i}.noDup.rl.bam"
	echo $cmd >> $out
	echo 'kill -USR2 %1' >> $out
	echo $cmd
	qsub -l mem=20G,time=40:: -o ${dir}/$log/calmd-${i}.o -e ${dir}/$log/calmd-${i}.e $out
	done
done

##index and sorting
#
samtool="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/samtools"
sampleList="sample.list"
dir=`pwd`
for a in `cat $sampleList`
do
	/bin/ls bam/$a/refine/*.noDup.rl.bam | cut -d "/" -f4 > $a.noDup.rl.bam.list
	echo "$a.noDup.rl.bam.list"
	for bam in `cat $a.noDup.rl.bam.list` 
	do
	i=`echo $bam |cut -d "." -f1 `
	echo $i
	log="bam/$a/refine/log"
	out=${dir}/$log/ixSrSt-${i}.sh
	echo '#!/bin/bash'  > $out
	echo 'uname -a' >> $out
	echo 'watch-myself &' >> $out
	cmd="$samtool sort ${dir}/bam/$a/refine/$bam ${dir}/bam/$a/refine/$i.noDup.rl.sorted"
	echo $cmd >> $out
	cmd="$samtool index ${dir}/bam/$a/refine/$i.noDup.rl.sorted.bam"
	echo $cmd >> $out
	cmd="$samtool idxstats ${dir}/bam/$a/refine/$i.noDup.rl.sorted.bam > ${dir}/bam/$a/refine/stat/$i.$a.noDup.rl.idxstats"
	echo $cmd >> $out
	cmd="$samtool flagstat ${dir}/bam/$a/refine/$i.noDup.rl.sorted.bam > ${dir}/bam/$a/refine/stat/$i.$a.noDup.rl.flagstat"
	echo $cmd >> $out
	echo 'kill -USR2 %1' >> $out
	echo $out
	qsub -l mem=6G,time=50:: -o ${dir}/$log/ixSrSt-${i}.o -e ${dir}/$log/ixSrSt-${i}.e $out
	done
done 

##fixmate and stat on mate
dir=`pwd`
setting=${dir}/global_setting_b37.sh
sampleList="sample.list.AC1"
. $setting
for a in `cat $sampleList`
do
	mkdir -p bam/$a
	mkdir -p bam/$a/refine
	mkdir -p bam/$a/refine/log
	for bam in `cat $a.noDup.rl.bam.sorted.list`
	do
	i=`echo $bam |cut -d "." -f1 `
	echo $i
	mkdir -p bam/$a/refine/stat/${i}
	log="bam/$a/refine/log"
	rm $log/fm-*
	out=${dir}/$log/fm-${i}.sh
	echo '#!/bin/bash'  > $out
	echo 'uname -a' >> $out
	echo 'watch-myself & watch_pid=$!' >> $out
	
	echo "java -Xmx5g -Djava.io.tmpdir=${myroot} -jar $FIXMATE INPUT=${myroot}/bam/$a/refine/$bam OUTPUT=${myroot}/bam/$a/refine/${i}.noDup.rl.sorted.fxmate.bam SO=coordinate VALIDATION_STRINGENCY=SILENT" >> $out
	cmd="java -Xmx5g -Djava.io.tmpdir=${myroot} -jar $FIXSTAT INPUT=${myroot}/bam/$a/refine/${i}.noDup.rl.sorted.fxmate.bam OUTPUT=${myroot}/bam/$a/refine/stat/${i}/${i}.noDup.rl.sorted.fxmate_details HISTOGRAM_FILE=${myroot}/bam/$a/refine/stat/${i}/${i}.noDup.rl.sorted.fxmate.hist.pdf REFERENCE_SEQUENCE=$REF"
	echo $cmd >> $out
		cmd="java -Xmx5g -Djava.io.tmpdir=${myroot} -jar $GCbias INPUT=${myroot}/bam/$a/refine/${i}.noDup.rl.sorted.bam OUTPUT=${myroot}/bam/$a/refine/stat/${i}/${i}.noDup.rl.sorted.fxmate.gc.bias_detail CHART=${myroot}/bam/$a/refine/stat/${i}/${i}.gcBias.pdf REFERENCE_SEQUENCE=$REF"
	echo $cmd >> $out
	cmd="samtools index ${dir}/bam/$a/refine/${i}.noDup.rl.sorted.fxmate.bam"
	echo $cmd >> $out
	#echo "gzip ${myroot}/bam/$a/refine/${i}.noDup.rl.bam" >> $out
	echo 'kill -USR2 $watch_pid; wait' >> $out
	echo $cmd
	qsub -l mem=6G,time=10:: -o ${dir}/$log/fm-${i}.o -e ${dir}/$log/fm-${i}.e $out 
	done
done

####################
###qc after reallining GARK
###################

dir=`pwd`
setting=${dir}/global_setting_b37.sh
sampleList="sample.list.AC1"
. $setting
for a in `cat $sampleList`
do
	/bin/ls bam/$a/refine/*.noDup.rl.sorted.bam | cut -d "/" -f4 > $a.noDup.rl.bam.sorted.list
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

##qc before
dir=`pwd`
. ${dir}/global_setting_b37.sh
sampleList="sample.list"
bamlist="bam.name.list"
Inputdir="/ifs/scratch/c2b2/ngs_lab/ngs/Projects/WG-seq/120309_ANDREA_GABRIELLE_6_HUMAN_GENOME_30X_PE_HISEQ/SAMPLES/AC1/BAM"

for a in `cat $sampleList`
do
	mkdir bam/$a
	mkdir bam/$a/
	mkdir bam/$a/log
	mkdir bam/$a/stat
	for bam in `cat $bamlist`
	do
	i=`echo $bam |cut -d "." -f1 `
	echo $i
	log="bam/$a/log"
	out=${dir}/$log/qc-${i}.sh
	echo '#!/bin/bash'  > $out
	echo 'uname -a' >> $out
	echo 'watch-myself & watch_pid=$!' >> $out
	cmd="java -Xmx4g -jar $GATKJAR -R $REF -knownSites $DBSNP -I ${Inputdir}/$bam -T CountCovariates -cov ReadGroupCovariate -cov QualityScoreCovariate -cov PositionCovariate -cov CycleCovariate -cov DinucCovariate -recalFile ${dir}/bam/$a/stat/$i.raw.qc.csv"
	echo $cmd >> $out
	cmd="java -Xmx4g -jar $GATKCOV -recalFile ${dir}/bam/$a/stat/$i.raw.qc.csv -outputDir ${dir}/bam/$a/stat/${i}/ -ignoreQ 5"
	echo $cmd >> $out
	echo 'kill -USR2 $watch_pid; wait' >> $out
	echo $cmd
	qsub -l mem=6G,time=10:: -o ${dir}/$log/qc-${i}.o -e ${dir}/$log/qc-${i}.e $out $setting
	done
done

###########DEATH of coverage
dir=`pwd`
setting=${dir}/global_setting_b37.sh
sampleList="sample.list.AC1"
. $setting
for a in `cat $sampleList`
do
	mkdir bam/$a
	mkdir bam/$a/refine
	mkdir bam/$a/refine/log
	for bam in `cat $a.noDup.rl.bam.sorted.list`
	do
	i=`echo $bam |cut -d "." -f1 `
	echo $i
	log="bam/$a/refine/log"
	out=${dir}/$log/dc-${i}.sh
	echo '#!/bin/bash'  > $out
	echo 'uname -a' >> $out
	echo 'watch-myself & watch_pid=$!' >> $out
	cmd="java -Xmx5g -jar $GATKJAR -R $REF -I ${dir}/bam/$a/refine/$bam  -T DepthOfCoverage -o ${dir}/bam/$a/refine/stat/$i.2.coverage -ct 1 -ct 5 -ct 10 -ct 15 -ct 20 -L ${i}"
	echo $cmd >> $out
	echo 'kill -USR2 $watch_pid; wait' >> $out
	echo $cmd
	qsub -l mem=6G,time=06:: -o ${dir}/$log/dc-${i}.o -e ${dir}/$log/dc-${i}.e $out
	done
done
#######################CALLING VARIANTS

##samtool callingVariants for sample
######################################################
dir=`pwd`
setting=global_setting_b37.sh
sampleList="sample.list.AC1"
. $setting
for a in `cat $sampleList`
do
	mkdir -p bam/$a/refine/calling
	mkdir -p bam/$a/refine/calling/log
	for bam in `cat $a.noDup.rl.bam.sorted.list`
	do
	i=`echo $bam |cut -d "." -f1 `
	echo $i
	log="bam/$a/refine/calling/log"
	out=${dir}/$log/calSNV-6-${i}.sh
	echo '#!/bin/bash'  > $out
	echo 'uname -a' >> $out
	echo 'watch-myself & watch_pid=$!' >> $out
	echo "source $myroot/$setting" >> $out
	cmd="samtools mpileup -Eug -d 400 -q 1 -C 50 -6 -f $REF ${dir}/bam/$a/refine/$bam | bcftools view -bvcg - > ${dir}/bam/$a/refine/calling/${i}.var.raw.6.bcf"
	echo $cmd >> $out
	cmd="bcftools view ${dir}/bam/$a/refine/calling/${i}.var.raw.6.bcf | vcfutils.pl varFilter -D 100 >  ${dir}/bam/$a/refine/calling/${i}.var.flt.6.vcf"
	echo $cmd >> $out
	echo 'kill -USR2 $watch_pid; wait' >> $out
	echo $cmd
	qsub -l mem=3G,time=5:: -o ${dir}/$log/calSNV-6-${i}.o -e ${dir}/$log/calSNV-6-${i}.e $out
	done
done



################################
##samtool callingVariants for paiered sample
###############################################
dir=`pwd`
setting=global_setting_b37.sh
sampleList="sample.list.all"
chrList="chrom.list"
. $setting
for a in `cat $sampleList`
do
	for b in `cat $sampleList`
	do
	if [ "$a" = "$b" ];
	then
	echo "same"
	else
	mkdir -p bam/callingJoin
	mkdir -p bam/callingJoin/log
	for bam in `cat AC1.noDup.rl.bam.sorted.list`
	do
	i=`echo $bam |cut -d "." -f1 `
	echo $i
	log="bam/callingJoin/log"
	out=${dir}/$log/calSNV-6.joint-${i}.sh
	echo '#!/bin/bash'  > $out
	echo 'uname -a' >> $out
	echo "source $myroot/$setting" >> $out
	echo 'watch-myself & watch_pid=$!' >> $out
	cmd="samtools mpileup -DSu -6 -f $REF ${dir}/bam/$a/refine/$bam ${dir}/bam/$b/refine/$bam | bcftools view -bvcgT pair - >  ${dir}/bam/callingJoin/${i}.joint.var.raw.6.bcf"
	echo $cmd >> $out
	cmd="bcftools view ${dir}/bam/callingJoin/${i}.joint.var.raw.6.bcf | vcfutils.pl varFilter -D 200 >  ${dir}/bam/callingJoin/${i}.joint.var.flt.6.vcf"
	echo $cmd >> $out
	echo 'kill -USR2 $watch_pid; wait' >> $out
	echo $cmd
	qsub -l mem=4G,time=15:: -o ${dir}/$log/calSNV.joint-6-${i}.o -e ${dir}/$log/calSNV.joint-6-${i}.e $out
	# jobid=`qsub ...    2>&1 | perl -ne '($jobid)=/(\d+)/; print $jobid'`
	done
	fi
	done
done

## merge vcf files) per samples 
dir=`pwd`
setting='global_setting_b37.sh'
sampleList="sample.list.all"
source ./$setting
for a in `cat $sampleList`
do
	mkdir -p bam/$a/refine/calling
	mkdir -p bam/$a/refine/calling/log
	mkdir -p bam/$a/refine/calling/final
	touch $myroot/$a.var.flt.6.vcf.gz.list
	/bin/ls bam/$a/refine/calling/*.var.flt.6.vcf > $a.var.flt.vcf.6.list
	for al in `cat $a.var.flt.vcf.6.list`
	do
	vcf=`basename $al`
	echo $vcf
	i=`echo $vcf |cut -d "." -f1 `
	echo $i
	cmd="bgzip $al"
	echo $cmd 
	$cmd
	cmd="tabix $al.gz"
	echo $cmd 
	$cmd
	/bin/ls bam/$a/refine/calling/*.var.flt.6.vcf.gz > $myroot/$a.var.flt.6.vcf.gz.list
	done
	cmd="vcf-concat -f $myroot/$a.var.flt.6.vcf.gz.list > $myroot/bam/$a/refine/calling/final/$a.flt.6.vcf.gz"
	echo $cmd 
	$cmd
	cmd="gunzip $myroot/bam/$a/refine/calling/final/$a.flt.6.vcf.gz"
	echo $cmd 
	$cmd
done

###annovar 
dir=`pwd`
setting='global_setting_b37.sh'
sampleList="sample.list.AC1"
source ./$setting
for a in `cat $sampleList`
do
	mkdir -p bam/$a/refine/calling/final/log
	log="bam/$a/refine/calling/final/log"
	out=${dir}/$log/an.sh
	echo '#!/bin/bash' >> $out
	echo 'uname -a' >> $out
	echo "source $myroot/$setting" >> $out
	echo 'watch-myself & watch_pid=$!' >> $out
	echo "convert2annovar.pl /${myroot}/bam/$a/refine/calling/final/$a.flt.vcf -format vcf4 -includeinfo -allallele >  ${myroot}/bam/$a/refine/calling/final/$a.flt.vcf.annovar" >> $out
	echo "summarize_annovar.pl  ${myroot}/bam/$a/refine/calling/final/$a.flt.vcf.annovar  /ifs/scratch/c2b2/ngs_lab/sz2317/softwares/annovar/humandb/  -outfile ${myroot}/bam/$a/refine/calling/final/$a.flt.vcf.annovar.summary -ver1000g 1000g2012feb -verdbsnp 132 -genetype=refgene --buildver hg19" >> $out
	echo 'kill -USR2 $watch_pid; wait' >> $out
	qsub -l mem=5G,time=50:: -o ${dir}/$log/an.o -e ${dir}/$log/an.e $out
done

## merge vcf files) joincalling
dir=`pwd`
setting='global_setting_b37.sh'
sampleList="sample.list.all"
source ./$setting
for a in `cat $sampleList`
do
	mkdir -p bam/callingJoin
	mkdir -p bam/callingJoin/log
	mkdir -p bam/callingJoin/final
	touch $myroot/callingJoin.var.flt.6.vcf.gz.list
	/bin/ls bam/callingJoin/*.var.flt.6.vcf > callingJoin.var.flt.vcf.6.list
	for al in `cat callingJoin.var.flt.vcf.6.list`
	do
	vcf=`basename $al`
	echo $vcf
	i=`echo $vcf |cut -d "." -f1 `
	echo $i
	cmd="bgzip $al"
	echo $cmd 
	$cmd
	cmd="tabix $al.gz"
	echo $cmd 
	$cmd
	/bin/ls bam/$a/refine/calling/*.var.flt.6.vcf.gz > $myroot/callingJoin.var.flt.6.vcf.gz.list
	done
	cmd="vcf-concat -f $myroot/callingJoin.var.flt.6.vcf.gz.list > $myroot/bam/$a/refine/calling/final/callingJoin.flt.6.vcf.gz"
	echo $cmd 
	$cmd
	cmd="gunzip $myroot/bam/$a/refine/calling/final/callingJoin.flt.6.vcf.gz"
	echo $cmd 
	$cmd
done


## downstream analysis after jointcalling
Include only sites with all Non-Reference Allele Counts within the specified range > 2

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

#######calling with varscan
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