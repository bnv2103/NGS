#!/bin/bash

fastadir=$1
outdir=$2

mkdir -p $outdir 

for g in `ls $fastadir`
do
    mkdir -p $outdir/$g
    for f in `ls $fastadir/$g/*fasta`
    do
        python /ifs/home/c2b2/ngs_lab/wsd2102/git_code/NGS/IonTorrent/covariate_table/fasta2homo.py $f $outdir/$g/`basename $f | cut -f1 -d "."`.homopickle
    done
done
