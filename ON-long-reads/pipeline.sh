#!/bin/sh
#$ -cwd

snprate=$1
errrate=$2
coverage=$3
maxreadlen=$4
minreadlen=$5
region=$6
iteration=$7
pval=$8

snprate=0.01
errrate=0.04
coverage=10
maxreadlen=5000
minreadlen=2000
region=$1
iteration=0
pval=1.1

base=$PWD
if [[ $region -eq 0 ]]
then
	pre=${snprate}_${errrate}_${coverage}_${maxreadlen}_${minreadlen}_${iteration}
else
	pre=${snprate}_${errrate}_${coverage}_${maxreadlen}_${minreadlen}_${iteration}
fi
read_base=${base}/reads/${pre}
vcf_base=${base}/vcfs/${pre}
output_base=${base}/HMM/output/${pre}

#gcc -g simulate_read.c -o simulate_read 

if [[ $region -eq 0 ]]
then
	let "qmem=400*$coverage"
	qsub -sync y -t 1-2 -l mem=4G,time=16:: ./simulate_read.sh $snprate $errrate $coverage $maxreadlen $minreadlen $region $iteration ${read_base}
	echo simulation complete

	let "qtime=3*$coverage"
	qsub -sync y -t 1-2 -l mem=8G,time=8:: ./long_read_map.sh ${read_base} $region
	rm ${read_base}.ibs
	for i in `seq 1 2`
	do
		cat ${read_base}_${i}.ibs >> ${read_base}.ibs
		#rm ${read_base}_${i}.ibs
	done
	echo alignment complete

	let "qtime=3*$coverage"
	qsub -sync y -t 1-2 -l mem=1G,time=8:: ./snp_calling.sh ${read_base} $pval ${vcf_base} $region
	echo snp calling complete

	cd HMM
	qsub -sync y -t 1-2 -l mem=32G,time=8:: ./program.sh ${read_base} ${vcf_base} $region $pval ${output_base}
	rm ${output_base}.${pval}.obs
	rm ${output_base}.${pval}.sta
	for i in `seq 1 2`
	do
		cat ${output_base}_${i}.${pval}.obs >> ${output_base}.${pval}.obs
		#rm ${output_base}_${i}.${pval}.obs
		cat ${output_base}_${i}.${pval}.sta >> ${output_base}.${pval}.sta
		#rm ${output_base}_${i}.${pval}.sta
	done
	echo HMM complete

	./analysis.sh ${output_base}.${pval} ${read_base}
else
	let "qmem=400*$coverage"
	qsub -sync y -l mem=4G,time=16:: ./simulate_read.sh $snprate $errrate $coverage $maxreadlen $minreadlen $region $iteration ${read_base}
	echo simulation complete

	let "qtime=3*$coverage"
	qsub -sync y -l mem=8G,time=8:: ./long_read_map.sh ${read_base} $region
	echo alignment complete

	let "qtime=3*$coverage"
	qsub -sync y -l mem=1G,time=8:: ./snp_calling.sh ${read_base} $pval ${vcf_base} $region
	echo snp calling complete

	cd HMM
	qsub -sync y -l mem=32G,time=8:: ./program.sh ${read_base} ${vcf_base} $region $pval ${output_base}
	echo HMM complete

	./analysis.sh ${output_base}_${region}.${pval} ${read_base}_${region}
fi

cd ..

echo Analysis complete

