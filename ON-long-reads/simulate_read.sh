#!/bin/sh
#$ -cwd

snprate=$1
errrate=$2
coverage=$3
maxreadlen=$4
minreadlen=$5
region=$6
iteration=$7
base=$8

if [[ $region -eq 0 ]]
then
	region=$SGE_TASK_ID
fi
prefix=${base}_${region}

echo -out_base=$prefix -ncells=$iteration -coverage=$coverage -snp_rate=$snprate -err_rate=$errrate -max_read_len=$maxreadlen -min_read_len=$minreadlen -region=$region

./simulate_read -out_base=$prefix -ncells=$iteration -coverage=$coverage -snp_rate=$snprate -err_rate=$errrate -max_read_len=$maxreadlen -min_read_len=$minreadlen -region=$region

# set args -out_base=0.01_0.04_10_5000_2000_1_0 -ncells=0 -coverage=10 -snp_rate=0.01 -err_rate=0.04 -max_read_len=5000 -min_read_len=2000 -region=1
# set args -out_base=/ifs/scratch/c2b2/ys_lab/aps2157/Haplotype/reads/0.01_0.04_10_5000_2000_0 -ncells=0 -coverage=10 -snp_rate=0.01 -err_rate=0.04 -max_read_len=5000 -min_read_len=2000 -region=0

# ./a.out -out_base=$prefix -ncells=$ind

