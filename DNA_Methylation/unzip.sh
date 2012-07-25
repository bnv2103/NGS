#!/bin/sh
#$ -cwd

for i in `ls *.gz`
do
	gunzip $i
	echo gunzip $i complete
done

