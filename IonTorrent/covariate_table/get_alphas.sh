#!/bin/bash

ctsdir=$1

for f in `ls $ctsdir/*.cts`
do
    mcode="path('`pwd`', path);
           path('/ifs/home/c2b2/ngs_lab/wsd2102/git_code/NGS/IonTorrent/lightspeed', path);
           path('/ifs/home/c2b2/ngs_lab/wsd2102/git_code/NGS/IonTorrent/fastfit', path);
           ctsfile = '`pwd`/$f';
           dat = load(ctsfile);
           alphai = 1000*[.002 .002 .002 .002];
           [~, fname, ~] = fileparts(ctsfile);
           alphai(find(fname(1) == 'ACGT')) = 0.994;
           alpha = polya_fit(dat, alphai);
           outfile = '`pwd`/$ctsdir/`basename $f .cts`.alpha';
           dlmwrite(outfile, alpha, '\t');
           exit;"
    echo matlab -nosplash -nodesktop -nodisplay -singleCompThread -nojvm -r \"$mcode\" | qsub -l mem=4G,time=1:: -o `pwd`/get_alphas.o -e `pwd`/get_alphas.e  
done
