#!/bin/sh
#$ -cwd
### 18G,1::

# Note: Please keep all fastq files in a single folder due to a restriction during the variant calling step.
# If you cannot do that, you will have to manually move the files after the merge step to a common folder and then do the variant calling

file=$1
file1=$2
file2=$3
chr=$4
build=$5

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
	qsub -l mem=32G,time=8:: ./run_gsnap.sh $i $chr
done

fi

while [[ -n `qstat | grep run_gsnap` ]]
do
	sleep 100
done

rm failed.list
for i in `ls run_gsnap.sh.e*`
do
	if [[ -z `grep Processed $i` ]]
	then
		head -4 $i | tail -1 | awk -F" " '{split($NF,file,".");print file[1];}' >> failed.list
	fi
done
if [[ -a failed.list ]]
then
	echo Rerun failed tests
	exit
fi

echo "GSNAP COMPLETE"

for i in `cat $file`
do
	goby 8g sort --output $i $i
done

echo "SORT COMPLETE"

for i in `seq 1 8`
do
	sample=`cat Sample_RK${i}.list`
	goby 16g concatenate-alignments $sample -o Sample_RK${i}
	goby 2g compact-file-stats Sample_RK${i}.entries
done

echo "MERGE COMPLETE"

merged_file="merged.list"

ls Sample_RK*.entries | cut -d'.' -f 1 > $merged_file

# An important restriction in this step is that the paths to the files should not contain folders. Use of folders requires the '/' character which is also being used in mentioning group names in the --compare switch
goby 30g discover-sequence-variants `cat ${merged_file}` --compare A/B --groups A=`cat ${file1} | tr ' ' ','`/B=`cat ${file2} | tr ' ' ','`  --format methylation --output meth.vcf --genome $chr

echo "METHYLATION COMPLETE"

