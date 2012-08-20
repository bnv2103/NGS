#!/bin/python

# count calls in mpileup and correlate with reference and homopolymer information
# for now it doesn't do anything with quality score


import pickle, sys, re, string
from numpy import mean

# mpileup file (for a specific amplicon)
fin = open(sys.argv[1], 'r')

# corresponding homo pickle file (for same amplicon)
homo = pickle.load(open(sys.argv[2], 'r'))

# what filtering do we want to do?? proabably depth? are we weighting observations by depth?
depth_thresh = 100
germline_thresh = .6 # up to 40% nonreference is considered not germline

for line in fin:
        # tokenize line in mpileup file
        [chr, pos, ref, depth, reads, Q, mQ] = line.split()

        # convert read string to simple sequence of nucleotides
        # capital letters for forward strand, lower case for reverse
        reads = re.sub(r'\^\^','', reads)
        reads = re.sub(r'\^.','', reads)
        reads = re.sub(r'\$','', reads)
        reads = re.sub('\.', ref, reads)
        reads = re.sub(',', string.lower(ref), reads)

        #remove indels (WAIT: we want to include this as covariate info right?)
        i = 0
        while i < len(reads):
                i = string.find(reads, '+', i, len(reads))
                if i < 0:
                        break
                match = re.search(r'\d+', reads[i+1:])
                insert_len = int(match.group(0))
                reads = reads[:i] + reads[i+insert_len+(match.end()-match.start()) + 1:]

        i = 0
        while i < len(reads):
                i = string.find(reads, '-', i, len(reads))
                if i < 0:
                        break
                match = re.search(r'\d+', reads[i+1:])
                delete_len = int(match.group(0))
                reads = reads[:i] + reads[i+delete_len+(match.end()-match.start()) + 1:]

        # make sure lengths are still consistent
        assert(len(Q) == len(reads) and int(depth) == len(reads) and len(mQ) == len(reads))

        count_table = {'forward':{'A':reads.count('A'), 'C':reads.count('C'), 'G':reads.count('G'), 'T':reads.count('T')}, \
                       'reverse':{'A':reads.count('a'), 'C':reads.count('c'), 'G':reads.count('g'), 'T':reads.count('t')}}

        depth_for = sum(count_table['forward'].itervalues())
        depth_rev = sum(count_table['reverse'].itervalues())
        if depth_for < depth_thresh or depth_rev < depth_thresh or \
           count_table['forward'][ref] < germline_thresh*depth_for or count_table['reverse'][ref] < germline_thresh*depth_rev:
            continue

        # generate seperate quality strings for each nucleotide and +/- 
        qA = ''.join([Q[i] for i, x in enumerate(reads) if x == "A"])
        qT = ''.join([Q[i] for i, x in enumerate(reads) if x == "T"])
        qC = ''.join([Q[i] for i, x in enumerate(reads) if x == "C"])
        qG = ''.join([Q[i] for i, x in enumerate(reads) if x == "G"])

        qa = ''.join([Q[i]  for i, x in enumerate(reads) if x == "a"])
        qt = ''.join([Q[i]  for i, x in enumerate(reads) if x == "t"])
        qc = ''.join([Q[i]  for i, x in enumerate(reads) if x == "c"])
        qg = ''.join([Q[i]  for i, x in enumerate(reads) if x == "g"])

        # generate seperate mapping quality strings for each nucleotide and +/- 
        mqA = ''.join([mQ[i] for i, x in enumerate(reads) if x == "A"])
        mqT = ''.join([mQ[i] for i, x in enumerate(reads) if x == "T"])
        mqC = ''.join([mQ[i] for i, x in enumerate(reads) if x == "C"])
        mqG = ''.join([mQ[i] for i, x in enumerate(reads) if x == "G"])

        mqa = ''.join([mQ[i]  for i, x in enumerate(reads) if x == "a"])
        mqt = ''.join([mQ[i]  for i, x in enumerate(reads) if x == "t"])
        mqc = ''.join([mQ[i]  for i, x in enumerate(reads) if x == "c"])
        mqg = ''.join([mQ[i]  for i, x in enumerate(reads) if x == "g"])

        # write data

        print chr+'\t'+pos+'\t'+ref+'\t'+ \
              str(count_table['forward']['A'])+'\t'+str(count_table['forward']['C'])+'\t'+str(count_table['forward']['G'])+'\t'+str(count_table['forward']['T'])+'\t'+ \
              qA+'\t'+qC+'\t'+qG+'\t'+qT+'\t'+\
              mqA+'\t'+mqC+'\t'+mqG+'\t'+mqT+'\t'+\
              str(homo[chr+':'+pos]['forward']) + '\t' + \
              str(count_table['reverse']['A'])+'\t'+str(count_table['reverse']['C'])+'\t'+str(count_table['reverse']['G'])+'\t'+str(count_table['reverse']['T'])+'\t'+ \
              qa+'\t'+qc+'\t'+qg+'\t'+qt+'\t'+\
              mqa+'\t'+mqc+'\t'+mqg+'\t'+mqt+'\t'+\
              str(homo[chr+':'+pos]['reverse'])

