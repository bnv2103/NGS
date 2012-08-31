#!/bin/python

# generates homopolymer pickle based on reference fasta (presumably from a single amplicon)
# this script generates the homopolymer data corresponding to a fasta (only one region ">" line)

import sys, pickle, re
from covariatefuncs import homodist

fasta = sys.argv[1]
homoout = sys.argv[2]

lines = open(fasta, 'r').readlines()
chr, start, end = re.split('[> : -]', lines[0][1:])
fullfasta = ''.join(lines[1:]).replace('\n', '')

outdict = {}
keys = [chr+':'+str(x) for x in range(int(start), int(end) + 1)]
for pos in range(len(fullfasta)):
    outdict[keys[pos]] = {'forward':homodist(fullfasta)[pos], 'reverse':homodist(fullfasta[::-1])[-pos-1]}

pickle.dump(outdict, open(homoout, 'w'))
