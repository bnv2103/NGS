#!/bin/bash
# -cwd

# this script runs the mpileupcts.py script for each amplicon in each gene in each sample 

# sample directory with mpileups
samples=`pwd`/$1
# homopickle directory
homodir=`pwd`/$2
# output directory (ctsandhomos)
datadir=`pwd`/$3
mkdir -p $datadir

for s in `ls -d $samples/*.mpileups`
do
    # sample directory in output
    sampledir=$datadir/`basename $s | cut -f1 -d "."`
    mkdir -p $sampledir 
    for g in `ls $s`
    do
        genedir=$datadir/`basename $s | cut -f1 -d "."`/$g
        mkdir -p $genedir 
        for f in `ls $s/$g/*mpileup`
        do
            echo python /ifs/home/c2b2/ngs_lab/wsd2102/git_code/NGS/IonTorrent/covariate_table/mpileupcts.py $f $homodir/$g/`basename $f | cut -f1 -d"."`.homopickle | qsub -l mem=1G,time=2:: -o $genedir/`basename $f | cut -f1 -d"."`.ctsandhomo -e `pwd`/get_cts.e
        done
    done
done
