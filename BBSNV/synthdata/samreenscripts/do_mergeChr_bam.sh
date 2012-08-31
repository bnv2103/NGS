#!/bin/bash
#$ -cwd


#Merge_sort_index the bam files which are present by chromosome names in all the directories which are found in the file list
# Input arg[0] - a file containing a list of directories which may contain bam files or symbolic links to bam files.  each directory should contain a header file created with samtools view -H, called "head.sam"
# Input arg[1] - output directory

# To Run: ./do_mergeChr_bam.sh AC234DirectoryList AC234

list=$1
first=`head -1 $list`

firstsg="$2/smfile"
samtools view -H $first > $firstsg

# firstsg=$first"/head.sam"

if [ ! -d logs ]; then mkdir logs; fi

grep "@HD" $firstsg > $2/outheader
grep "@SQ" $firstsg >> $2/outheader
if [ -e $2/temp ] ; then rm  $2/temp;  fi
for i in `cat $list` ;do
	samtools view -H $i | grep "@RG" >>  $2/temp
done
sort -u  $2/temp >> $2/outheader
grep "@PG" $firstsg >> $2/outheader

rm $2/temp

for i in `ls $first/*.bam ` ;do
	chr=`basename $i`
	CMD="qsub -o `pwd`/$chr.o -e `pwd`/$chr.e -N merge.$chr -l mem=10G,time=15:: $WGS/mergeChr_bam.sh $chr `readlink -f $list` $2/outheader `readlink -f $2`/"
	echo $CMD
	$CMD
done
