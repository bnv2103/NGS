#!/bin/python

# generate cell data files

import sys, os

infile = open(sys.argv[1], 'r')
outdir = sys.argv[2]
os.mkdir(outdir)

comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}

for line in infile:

    chr, pos, ref, A, C, G, T, qA, qC, qG, qT, mqA, mqC, mqG, mqT, H, a, c, g, t, qa, qc, qg, qt, mqa, mqc, mqg, mqt, h = line.rstrip('\n').split('\t')

    outfile_for = open(outdir+'/'+ref+'_'+H.lstrip('[').rstrip(']').replace(', ', '_')+'.cts', 'a')
    outfile_rev = open(outdir+'/'+comp[ref]+'_'+h.lstrip('[').rstrip(']').replace(', ', '_')+'.cts', 'a')

    outfile_for.write('\t'.join([A, C, G, T])+'\n')
    outfile_rev.write('\t'.join([t, g, c, a])+'\n')

    outfile_for.close()
    outfile_rev.close()
