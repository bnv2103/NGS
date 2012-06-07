#!/bin/bash
#$ -cwd
# inputs: fastq file directory

# file1=$1
# file2=$1
 dir="./"

# cd $dir 

# echo $file1
# echo $file2
# echo $dir

for f in *120418*fastq
do  
    
    for g in *120509*fastq
    do
	# echo $f
	ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/checkName.rb $f $g $dir

    done

done