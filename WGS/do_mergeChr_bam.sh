#!/bin/bash
#$ -cwd


#Merge_sort_index the bam files which are present by chromosome names in all the directories which are found in the file list
# Input arg[0] - a file containing a list of directories which may contain bam files or symbolic links to bam files.  each directory should contain a header file created with samtools view -H, called "head.sam"
# Input arg[1] - output directory

# To Run: ./do_mergeChr_bam.sh AC234DirectoryList AC234

list=$1
first=`head -1 $list`
firstsg=$first"/head.sam"

if [ ! -d logs ]; then mkdir logs; fi

grep "@HD" $firstsg > $2/outheader
grep "@SQ" $firstsg >> $2/outheader
if [ -e temp ] ; then rm temp; fi
for i in `cat $list` ;do
	grep "@RG" $i"/head.sam" >> temp
done
sort -u temp >> $2/outheader
grep "@PG" $firstsg >> $2/outheader


for i in `ls $first/*.bam ` ;do
	chr=`basename $i`
    CMD="qsub -o logs/$chr.o -e logs/$chr.e -N merge.$chr -l mem=5G,time=6:: $WGS/mergeChr_bam.sh $chr `readlink -f $list` outheader `readlink -f $2`/"
	echo $CMD
	$CMD
    exit
done
