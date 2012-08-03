#!/bin/sh
#$ -cwd
### mem=2G,time=2::

base=$1

sample=${base}${SGE_TASK_ID}

goby 1g alignment-to-counts $sample -o $sample

