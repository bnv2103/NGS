#!/bin/bash
#$ -cwd
#$ -l mem=2G,time=2::

SAMTOOL="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/samtools"
list_bam=$1  ## list of deduped bams

for i in `cat list_bam`;
do
#	j=`basename $i`
	$SAMTOOL flagstat $i > $i.flagstat
done

# coverage= sum of all mapped reads /size of human genome

