#!/bin/bash
#$ -cwd
file1=$1
file2=$2
dir=$3

cd $dir 

for f in $file1
do  
    
    for g in $file2
    do
       
	ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/checkName.rb $f $g $dir

    done

done