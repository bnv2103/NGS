#!/bin/python

# generate cell data files

import sys, os

# input file is the output file from get_cts.sh
infile = open(sys.argv[1], 'r')
# specify new output directory for covariate cell data
outdir = sys.argv[2]
os.mkdir(outdir)

# reverse complement hash
comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

# loop through all sites in the input file and write the counts vector to the appropriate covariate cell file (determined by reference base and homopolymer info)
for line in infile:

    chr, pos, ref, A, C, G, T, qA, qC, qG, qT, mqA, mqC, mqG, mqT, H, a, c, g, t, qa, qc, qg, qt, mqa, mqc, mqg, mqt, h = line.rstrip('\n').split('\t')

    # covariate cell file for forward reads, determined by reference base and forward homopolymer data
    outfile_for = open(outdir+'/'+ref+'_'+H.lstrip('[').rstrip(']').replace(', ', '_')+'.cts', 'a')
    # covariate cell file for reverse reads, determined by complement reference base and reverse homopolymer data
    outfile_rev = open(outdir+'/'+comp[ref]+'_'+h.lstrip('[').rstrip(']').replace(', ', '_')+'.cts', 'a')

    # append count vector for forward reads to forward read covariate cell file
    outfile_for.write('\t'.join([A, C, G, T])+'\n')
    # append complemented count vector for reverse reads to reverse read covariate cell file 
    outfile_rev.write('\t'.join([t, g, c, a])+'\n')

    # close the files
    outfile_for.close()
    outfile_rev.close()
