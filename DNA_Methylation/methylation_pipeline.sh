#!/bin/sh
#$ -cwd
### 18G,1::

file=$1
file1=$2
file2=$3
chr=$4
build=$4

if [[ $build -eq 1 ]]
then
	gmap_build -d $chr -k 15 --basesize=15 ${chr}/${chr}.fasta
	echo "GMAP BUILD COMPLETE"

	goby 5g build-sequence-cache --basename $chr ${chr}/${chr}.fasta
	echo "BUILD SEQUENCE CACHE COMPLETE"

	cmetindex -d $chr -k 15 -b 15
	echo "CMETINDEX COMPLETE"


for i in `cat $file`
do
	qsub -l mem=32G,time=2:: ./run_gsnap.sh $i $chr
done

while [[ -n `qstat | grep run_gsnap` ]]
do
	sleep 100
done

echo "GSNAP COMPLETE"

for i in `cat $file`
do
	goby 1g sort --output $i $i
done

echo "SORT COMPLETE"

for i in `seq 1 8`
do
	sample=`cat Sample_RK${i}.list`
	goby 4g concatenate-alignments $sample -o Sample_RK${i}
	goby 1g compact-file-stats Sample_RK${i}.entries
done

echo "MERGE COMPLETE"
fi

merged_file="merged.list"

ls Sample_RK*.entries | cut -d'.' -f 1 > $merged_file

goby 30g discover-sequence-variants `cat ${merged_file}` --compare A/B --groups A=`cat ${file1} | tr '\n' ','`/B=`cat ${file2} | tr '\n' ','`  --format methylation --output meth.vcf --genome $chr

echo "METHYLATION COMPLETE"

