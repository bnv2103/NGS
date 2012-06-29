#!/bin/python

# filters mpileup file based on % top nonref allele and min number of top nonref on forward and reverse
# annotates with quality, mapping quality, and homopolymer distance and size, and indel proximity, all segregated by base.

import sys, re, string
from numpy import mean, std
from scipy import stats

fin = open(sys.argv[1], 'r')
#fout = open(sys.argv[1] + '.can', 'w')

# filter thresholds on top nonreference allele frequency and count
freq_thresh = 0.01
count_thresh = 10

# filter for p-value in U test, for and rev
# for now we are keeping them all
# if both are low, that indicated false positive
p_thresh = 0#0.01

for line in fin:
        # tokenize line in mpileup file
        [chr, pos, ref, depth, reads, Q] = line.split()

        # convert read string to simple sequence of nucleotides
        # capital letters for forward strand, lower case for reverse
        reads = re.sub(r'\^\^','', reads)
        reads = re.sub(r'\^.','', reads)
        reads = re.sub(r'\$','', reads)
        reads = re.sub('\.', ref, reads)
        reads = re.sub(',', string.lower(ref), reads)

        #remove indels (WAIT: we want to include this right?)
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
        assert(len(Q) == len(reads) and int(depth) == len(reads))

        count_table = {'forward':{'A':reads.count('A'), 'C':reads.count('C'), 'G':reads.count('G'), 'T':reads.count('T')}, \
                       'reverse':{'A':reads.count('a'), 'C':reads.count('c'), 'G':reads.count('g'), 'T':reads.count('t')}}

        depth_for = sum(count_table['forward'].itervalues())
        depth_rev = sum(count_table['reverse'].itervalues())
        
        topnonrefct_for = max([n for b, n in count_table['forward'].iteritems() if b != ref])
        topnonrefct_rev = max([n for b, n in count_table['reverse'].iteritems() if b != ref])

        topnonref = ''
        for n, ct in count_table['forward'].iteritems():
            if ct == topnonrefct_for and n != ref:
                topnonref += n
        topnonref_rev = ''
        for n, ct in count_table['reverse'].iteritems():
            if ct == topnonrefct_rev and n != ref:
                topnonref_rev += n

        # top non reference nucleotides must agree between forward and reverse, also filter on count
        if topnonref != topnonref_rev or topnonrefct_for < freq_thresh*depth_for or topnonrefct_rev < freq_thresh*depth_rev or topnonrefct_for < count_thresh or topnonrefct_rev < count_thresh:
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

        q_table = {'forward':{'A':map(ord, qA), 'C':map(ord, qC), 'G':map(ord, qG), 'T':map(ord, qT)}, \
                   'reverse':{'A':map(ord, qa), 'C':map(ord, qc), 'G':map(ord, qg), 'T':map(ord, qt)}}



        # missing functionality:
        # get homopolymer info
        # get mapping quality
        # get indel info

        # NOTE: these require dealing with sam file, entire reads

        # write data
        #fout.write(chr+'\t'+pos+'\t'+ref+'\t'+qA+'\t'+qC+'\t'+qG+'\t'+qT+'\t'+qa+'\t'+qc+'\t'+qg+'\t'+qt+'\n')
        # fout.write(chr+'\t'+pos+'\t'+ref+'\t'+topnonref+'\t'+ \
        #           str(count_table['forward'][ref])+'\t'+str(count_table['forward'][topnonref])+'\t'+ \
        #           str(mean(q_table['forward'][ref]))+'\t'+str(mean(q_table['forward'][topnonref]))+'\t'+\
        #           str(count_table['reverse'][ref])+'\t'+str(count_table['reverse'][topnonref])+'\t'+ \
        #           str(mean(q_table['reverse'][ref]))+'\t'+str(mean(q_table['reverse'][topnonref]))+'\n')


        # using mann-whitney U

        u_for, p_for = stats.mannwhitneyu(q_table['forward'][ref], q_table['forward'][topnonref], use_continuity=False)
        u_rev, p_rev = stats.mannwhitneyu(q_table['reverse'][ref], q_table['reverse'][topnonref], use_continuity=False)

        # null hypothesis corresponds to mutation
        # mult p-val by 2 to get two-sided p-value
        p_for = 2*p_for
        p_rev = 2*p_rev

        if p_for < p_thresh or p_rev < p_thresh: continue

        print chr+'\t'+pos+'\t'+ref+'\t'+topnonref+'\t'+ \
              str(count_table['forward']['A'])+'\t'+str(count_table['forward']['C'])+'\t'+str(count_table['forward']['G'])+'\t'+str(count_table['forward']['T'])+'\t'+ \
              str(mean(q_table['forward'][ref]))+'\t'+ \
              str(std(q_table['forward'][ref]))+'\t'+ \
              str(mean(q_table['forward'][ref]))+'\t'+ \
              str(std(q_table['forward'][ref]))+'\t'+ \
              str(p_for)+'\t'+\
              str(count_table['reverse']['A'])+'\t'+str(count_table['reverse']['C'])+'\t'+str(count_table['reverse']['G'])+'\t'+str(count_table['reverse']['T'])+'\t'+ \
              str(mean(q_table['reverse'][ref]))+'\t'+ \
              str(std(q_table['reverse'][ref]))+'\t'+ \
              str(mean(q_table['reverse'][ref]))+'\t'+ \
              str(std(q_table['reverse'][ref]))+'\t'+ \
              str(p_rev)
