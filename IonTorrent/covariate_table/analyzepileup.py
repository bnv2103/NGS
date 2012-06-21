#!/bin/python

from re import compile, finditer, findall, match
from Bio.Seq import Seq
from covariatefuncs import homodist, cigar2align
import pickle, sys, os
import scipy.io

# full bam
bamname = sys.argv[1]
# ref
refname = sys.argv[2]
# region sam file, for amplicon
sam = open(sys.argv[3], 'r') 
samlines = sam.read().splitlines()
sam.close()
# bed file defining region
bed = open(sys.argv[4], 'r')
chr, start, end = bed.readline().split()
bed.close()

# the covariate hash table
table = pickle.load(open(sys.argv[5], 'r'))

# kludge to use matlab's cigar2align
alignedf = open(sys.argv[3]+'.alignedMATLAB', 'r')
aligned_data = alignedf.read().splitlines()
alignedf.close()

outdir = sys.argv[4]+'.dataforMLE'
os.mkdir(outdir)

# looping over each pileup
for i in range(int(start), int(end) + 1):
    # make temp sam file of reads covering specific position i
    f = os.popen('samtools view '+bamname+' '+chr+':'+str(i)+'-'+str(i))
    # these hold counts in the covariate table for the pileup bases
    tablecounts_for = []
    tablecounts_rev = []
    # these hold total pileup base counts
    counts_for = {'A':0, 'C':0, 'G':0, 'T':0}
    counts_rev = {'A':0, 'C':0, 'G':0, 'T':0}
    # looping over each pileup read
    true = ''
    for line in f:
        qname, flag, rname, pos = line.rstrip('\n').split()[:4]
        flag = int(flag)
        reverse_flag = bool(flag & 16)
        pos = int(pos)
        for samlinenum in range(len(samlines)):
            if samlines[samlinenum] == line.rstrip('\n'):
                regionindex = samlinenum
                break
        aligned_read, aligned_quality = aligned_data[regionindex].split() 

        f2 = os.popen('samtools faidx '+refname+' '+rname+':'+str(pos)+'-'+str(pos + len(aligned_read) - 1))
        aligned_ref = ''.join([line2.rstrip('\n') for line2 in f2][1:])
        f2.close()

        positions = range(pos, pos + len(aligned_ref))

        # reverse complement
        if reverse_flag:
            aligned_read = str(Seq(aligned_read).reverse_complement());
            aligned_ref = str(Seq(aligned_ref).reverse_complement());
            aligned_quality = aligned_quality[::-1];
            positions = positions[::-1]

        hs_all, hd_all = homodist(aligned_ref)

        # index into the read to get this position
        siteindex = positions.index(i)

        ref = aligned_ref[siteindex]        
        if not true and not reverse_flag:
            true = ref
        call = aligned_read[siteindex]
        hs = hs_all[siteindex]
        hd = hd_all[siteindex]
        q = aligned_quality[siteindex]
        # q = ord(aligned_quality[siteindex]) - 33

        if call == '-':
            continue

        # lookup counts in the covariate table
        # NOTE: this assumes the cell exists!
        new_tablecount = table['\t'.join([ref, q, str(hs), str(hd)])]
        # new_tablecount_total = sum(table_counts.itervalues())

        ## for now we skip bases in the pileup that don't have large occupancy in covariate table
        #if table_counts_total < 1000:
        #    continue

        if reverse_flag:
            counts_rev[call] += 1
            tablecounts_rev.append([new_tablecount['A'], new_tablecount['C'], new_tablecount['G'], new_tablecount['T']])
        else:
            counts_for[call] += 1
            tablecounts_for.append([new_tablecount['A'], new_tablecount['C'], new_tablecount['G'], new_tablecount['T']])

    f.close()

    # ok, now we have the stack of p_hat and base counts for forward and reverse seperately, time to do dirichlet MLE 
    # for now the easiest approach is to use Minka's matlab code, but should be easily portable to python numpy 

    # below we are writing a .mat file for each pileup

    outdict = {'counts_for':[counts_for['A'], counts_for['C'], counts_for['G'], counts_for['T']], \
               'counts_rev':[counts_rev['A'], counts_rev['C'], counts_rev['G'], counts_rev['T']], \
               'tablecounts_for':tablecounts_for, \
               'tablecounts_rev':tablecounts_rev}

    scipy.io.savemat(outdir+'/'+chr+'_'+str(i)+'_'+true+str(Seq(true).reverse_complement()), outdict)

