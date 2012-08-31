#!/bin/python

from re import compile, finditer, findall, match
from Bio.Seq import Seq
from covariatefuncs import homodist, cigar2align
import sys, os

sam = open(sys.argv[1], 'r') 
out = open(sys.argv[1]+'.errorpatterns', 'w')

# kludge to use matlab's cigar2align
alignedf = open('sample.sam.alignedMATLAB', 'r')

for line in sam:
    qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual = line.split()[:11]

    flag = int(flag)
    pos = int(pos)
    mapq = int(mapq)
    pnext = int(pnext)
    tlen = int(tlen)

    #aligned_read = alignedf.readline().rstrip('\n')

    if (flag & 4) or not match('chr[1-9 10-22 X Y]', rname) or pos < 1 or seq == '*' or qual == '*':
        continue

    reverse_flag = bool(flag & 16)
                  

    # this shows deletions as hyphens and omits insertions. Later if we want

    # using matlab kludge
    # aligned_read, aligned_quality = cigar2align(cigar, seq, qual)
    aligned_read, aligned_quality = alignedf.readline().split()
 
    pos = [pos + idx for idx in range(len(aligned_read))]
    # this is where I need the subprocess check_output to call 'samtools faidx sys.argv[2] '+rname+':'+str(pos[0])+'-'+str(pos + len(aligned_read) - 1)
    f = os.popen('samtools faidx '+sys.argv[2]+' '+rname+':'+str(pos[0])+'-'+str(pos[-1]))
    aligned_ref = ''.join([line.rstrip('\n') for line in f][1:])
    f.close()

    # reverse complement
    if reverse_flag:
        aligned_read = Seq(aligned_read).reverse_complement();
        aligned_ref = Seq(aligned_ref).reverse_complement();
        aligned_quality = aligned_quality[::-1];

    # ignore deletions from ref (insertions were removed by cigar2align)
    # dels = aligned_read ~= '-';
    # aligned_read = aligned_read(dels);
    # aligned_ref = aligned_ref(dels);
    # aligned_quality = aligned_quality(dels);
    # pos = pos(dels);
    
    # compute the error pattern variables

    GC = (aligned_ref.count('G') + aligned_ref.count('C'))/float(len(aligned_ref))
    hs, hd = homodist(str(aligned_ref))
    readlen = len(seq)
    for i in range(len(aligned_read)):
        if aligned_read[i] == '-': continue
        #region = rname[3:]+':'+str(pos[i])+'-'+str(pos[i])
        #f = os.popen('/ifs/scratch/c2b2/ngs_lab/sz2317/bin/src/tabix-0.2.5/tabix /ifs/scratch/c2b2/ngs_lab/ngs/Projects/IonTorrent/will/xunhai/covariate_table/dbSNP132b37.vcf.gz '+region)
        #if f.readline() == '':
        out.write('\t'.join(map(str,[rname, pos[i], aligned_ref[i], aligned_read[i], readlen, int(reverse_flag), aligned_quality[i], pos[i], hs[i], hd[i], GC]))+'\n')
        #f.close()

out.close()
