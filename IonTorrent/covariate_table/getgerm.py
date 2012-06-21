#!/bin/python

import sys
from re import match
fin = open(sys.argv[1], 'r')
fout = open(sys.argv[1]+'GERM', 'w')
for l in fin:
    chro, pos, ref, depth, pile, qual = l.split()
    if not match('chr[1-9 10-22 X Y]', chro) or not match('[A T C G]', ref): continue
    if (l.count(',') + l.count('.'))/float(depth) < .40:
        fout.write(chro+'\t'+pos+'\n')
