#!/bin/sh
#$ -cwd

# One-time calling
#bwa index -a bwtsw bcm_hg18.fasta

# Defaults are: -b 3 -q 5 -r 2
# Expt range: -b [4-6] -q [2-4] -r [1-2]

bwa bwasw -b 4 -q 4 -r 2 bcm_hg18.fasta simulated_reads.fq > simulated_chr20_b4q4r2.sam

