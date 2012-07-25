#!/bin/sh
#$ -cwd
## 1G,8::

prefix=$1
pval=$2
suffix=$3
region=$4

if [[ $region -eq 0 ]]
then
	region=$SGE_TASK_ID
fi
rbase=${prefix}_$region
vbase=${suffix}_$region

qsub -sync y -l mem=1G,time=8:: -t 1-10 mpileup.sh ${rbase} $region $vbase $pval

grep ^# ${vbase}.1.temp.vcf > ${vbase}.${pval}.vcf
mkdir -p ./tmp_${region}

grep -hv ^# ${vbase}.1.temp.vcf ${vbase}.2.temp.vcf ${vbase}.3.temp.vcf ${vbase}.4.temp.vcf ${vbase}.5.temp.vcf ${vbase}.6.temp.vcf ${vbase}.7.temp.vcf ${vbase}.8.temp.vcf ${vbase}.9.temp.vcf ${vbase}.10.temp.vcf | sort -u -k2,2g -T ./tmp_${region} | awk -F'\t' '{if($4!="N"&&length($5)==1) { split($8,info,";");split(info[1],dp,"="); if(dp[2]>=2) print $0;} }' >> ${vbase}.${pval}.vcf

rm -rf ./tmp_${region}
rm ${vbase}.bcf.1 ${vbase}.bcf.2 ${vbase}.bcf.3 ${vbase}.bcf.4 ${vbase}.bcf.5 ${vbase}.bcf.6 ${vbase}.bcf.7 ${vbase}.bcf.8 ${vbase}.bcf.9 ${vbase}.bcf.10
rm ${vbase}.1.temp.vcf ${vbase}.2.temp.vcf ${vbase}.3.temp.vcf ${vbase}.4.temp.vcf ${vbase}.5.temp.vcf ${vbase}.6.temp.vcf ${vbase}.7.temp.vcf ${vbase}.8.temp.vcf ${vbase}.9.temp.vcf ${vbase}.10.temp.vcf

