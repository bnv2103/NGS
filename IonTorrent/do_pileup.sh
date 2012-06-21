#!/bin/bash
#$ -cwd

bam=$1
#out=$2

ref="/ifs/home/c2b2/ngs_lab/ngs/data/resources/ionTorrent_hg19/hg19.fasta"

region="chr3:178935893-178936293"


samtools mpileup -BD -d 500 -m 1 -F 0.00001  -f $ref -r $region $bam # > $out
