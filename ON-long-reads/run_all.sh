#!/bin/sh
#$ -cwd
### mem=6G,time=16::

prefix=$1
ind=$2

gcc -g simulate_read.c
qsub -sync y -l mem=1G,time=1:: ./simulate_read.sh $prefix $ind
echo simulation complete
qsub -sync y -l mem=6G,time=1:: ./long_read_map.sh ${prefix}_${ind}
echo alignment complete
qsub -sync y -l mem=1G,time=1:: ./snp_calling.sh ${prefix}_${ind}
echo snp calling complete
