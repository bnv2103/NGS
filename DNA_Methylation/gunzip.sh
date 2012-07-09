#!/bin/sh
#$ -cwd

for i in `seq 1 8`
do
	cd Sample_RK$i
	num=`ls -lrt *.gz | wc -l`
	let "time=num/5"
	cp ../unzip.sh .
	qsub -l mem=2G,time=${time}:: unzip.sh
	cd ..
done

