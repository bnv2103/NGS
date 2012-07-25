#!/bin/bash
# -cwd

ctsandhomodir=$1
celldir=$2
outdir=$3

mkdir -p $outdir

for s in `ls  $ctsandhomodir`
do
    mkdir -p $outdir/$s 
    for g in `ls $ctsandhomodir/$s`
    do
        mkdir -p $outdir/$s/$g
        for f in `ls $ctsandhomodir/$s/$g/*ctsandhomo`
        do
            mcode="path('/ifs/home/c2b2/ngs_lab/wsd2102/git_code/NGS/IonTorrent/covariate_table', path);
                   path('/ifs/home/c2b2/ngs_lab/wsd2102/git_code/NGS/IonTorrent/lightspeed', path);
                   path('/ifs/home/c2b2/ngs_lab/wsd2102/git_code/NGS/IonTorrent/fastfit', path);
                   ctsandhomo2pvals('`pwd`/$f', '`pwd`/$celldir', '`pwd`/$outdir/$s/$g/`basename $f .ctsandhomo`.pvals');
                   exit"
            echo matlab -singleCompThread -nojvm -nodisplay -nosplash -r \"$mcode\" | qsub -l mem=2G,time=2:: -o `pwd`/get_pvals.o -e `pwd`/get_pvals.e
        done
    done
done
