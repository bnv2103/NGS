#!/bin/bash
#$ -cwd
#$ -l mem=1G,time=2::

inp=$1
SAMTOOLS="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/samtools"

if [ ! -d $inp"_split" ]; then mkdir $inp"_split"; fi

$SAMTOOLS view -H $inp > $inp"_split/"head.sam

#example of head.sam
#@HD     VN:1.0  GO:none SO:coordinate
#@SQ     SN:1    LN:249250621
#@SQ     SN:2    LN:243199373
#@SQ     SN:3    LN:198022430

while read line
do
 if [[ $line == "@SQ"* ]];then
	chr=`echo $line |cut -f2 -d" " | sed 's/^SN://'`
	echo "Extracting chromosome $chr"
	$SAMTOOLS view -hb $inp $chr >  $inp"_split/"$chr".bam" 
 fi
done <  $inp"_split/"head.sam

