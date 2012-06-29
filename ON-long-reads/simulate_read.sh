#!/bin/sh
#$ -cwd
### mem=1G,1::

prefix=$1
ind=$2
cov=$3

# ./simulate_read -out_base=$prefix -ncells=$ind -coverage=$cov
 ./a.out -out_base=$prefix -ncells=$ind

