#!/bin/python

# this script generates the homopolymer data corresponding to a fasta made by get_fasta.sh (only one region ">" line)

import sys, pickle, re
from covariatefuncs import homodist

fasta = sys.argv[1]
homoout = sys.argv[2]

lines = open(fasta, 'r').readlines()
chr, start, end = re.split('[> : -]', lines[0][1:])
homodat = homodist(''.join(lines[1:]).replace('\n', ''))

outdict = {}
keys = [chr+':'+str(x) for x in range(int(start), int(end) + 1)]
for pos in range(len(homodat['forward'])):
    outdict[keys[pos]] = {'forward':homodat['forward'][pos], 'reverse':homodat['reverse'][pos]}

pickle.dump(outdict, open(homoout, 'w'))
