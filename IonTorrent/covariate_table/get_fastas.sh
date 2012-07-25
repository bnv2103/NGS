#!/bin/bash
#$ -cwd

fasta=$1 # reference file
bed=$2 # a bed file listing genes and amplicon regions
outdir=$3
rm -rf $outdir
mkdir -p $outdir

N=`grep chr $bed | wc -l`
n=0

echo "generating fasta for each amplicon"

while read line; do
	if [ ${line:0:3} != chr ]; then continue; fi
	line=($line)
	chrom=${line[0]}
	begin=${line[1]}
	end=${line[2]}

	region=$chrom:$begin-$end
	amplicon=${line[3]}
	gene=${line[5]}

	genedir=$outdir/$gene
	mkdir -p $genedir
	outfile=$genedir/$amplicon.fasta

        samtools faidx $fasta $region >> $outfile

	let "n++"
	echo "$n of $N amplicons completed"
done < $bed
