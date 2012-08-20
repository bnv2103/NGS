#!/bin/bash
# -cwd


datadir=`pwd`/$1
mkdir -p $datadir

for s in `ls -d $datadir/../../sample_files/*.mpileups`
do
    sampledir=$datadir/`basename $s | cut -f1 -d "."`
    mkdir -p $sampledir 
    for g in `ls $s`
    do
        genedir=$datadir/`basename $s | cut -f1 -d "."`/$g
        mkdir -p $genedir 
        for f in `ls $s/$g/*mpileup`
        do
            echo python /ifs/home/c2b2/ngs_lab/wsd2102/git_code/NGS/IonTorrent/covariate_table/mpileupcts.py $f $datadir/../homopickles/$g/`basename $f | cut -f1 -d"."`.homopickle | qsub -l mem=1G,time=2:: -o $genedir/`basename $f | cut -f1 -d"."`.ctsandhomo -e `pwd`/get_ctsandhomo.e
        done
    done
done
