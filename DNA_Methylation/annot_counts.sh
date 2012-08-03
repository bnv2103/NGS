#!/bin/sh
#$ -cwd
### mem=2G,time=2::

# One time calling only

goby 1g annotations-to-counts mm10-capture-regions.tsv -o capture-regions

