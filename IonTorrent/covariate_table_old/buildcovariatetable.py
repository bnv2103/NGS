#!/bin/bash

import sys
import pickle

# load errorpattern file with germline sites excluded

table = {}
for line in open(sys.argv[1], 'r'):
    datline = line.rstrip('\n').split()
    # only looking at q, hs, hd
    cell = '\t'.join([datline[2], datline[6], datline[8], datline[9]])
    call = datline[3] 
    if cell not in table:
        table[cell] = {'A':0, 'C':0, 'G':0, 'T':0}
    table[cell][call] += 1 

pickle.dump(table, open(sys.argv[1]+'.covariatetable', 'w'))

