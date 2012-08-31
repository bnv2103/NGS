#!/bin/bash
#$ -cwd

bam=$1 # alignment file
bed=$2 # a bed file listing genes and amplicon regions
ref=$3 #"/ifs/home/c2b2/ngs_lab/ngs/data/resources/ionTorrent_hg19/hg19.fasta" # reference genome

outdir=$bam.mpileups
rm -rf $outdir
mkdir -p $outdir

N=`grep chr $bed | wc -l`
n=0

# threshold for mapping quality (note, we keep all base qualities)
q=20 #30

echo "generating pileups, mapQ filter = $q"

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
	outfile=$genedir/$amplicon.mpileup

	# samtools mpileup -BD -d 500 -m 1 -F 0.00001  -f $ref -r $region $bam >> $outfile
        samtools mpileup -s -Q 0 -q $q -f $ref -r $region $bam >> $outfile

	let "n++"
	echo "$n of $N amplicons completed"
done < $bed
