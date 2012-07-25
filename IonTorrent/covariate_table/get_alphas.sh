#!/bin/bash

ctsdir=$1

for f in `ls $ctsdir/*.cts`
do
    command="\"path('`pwd`', path); get_alphas('`pwd`/$f', '`pwd`/$ctsdir/`basename $f .cts`.alpha'); exit;\""
    echo matlab -nosplash -nodesktop -nodisplay -singleCompThread -nojvm -r "$command" | qsub -l mem=1G,time=1:: -o `pwd`/get_alphas.o -e `pwd`/get_alphas.e  
done
