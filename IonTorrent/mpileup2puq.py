#!/bin/python

# generates pileup quality file (puq) from mpileup input
import sys, re, string

fin = open(sys.argv[1], 'r')
fout = open(sys.argv[1] + '.puq', 'w')

#fout.write('## reference\tforwardA\tforwardC\tforwardG\tforwardT\treverseA\treverseC\treverseG\treverseT\n')
for line in fin:
	words = line.split()
	Q = words[5]
	reads = words[4] 
	reads = re.sub(r'\^\^','', reads)
        reads = re.sub(r'\^.','', reads)
        reads = re.sub(r'\$','', reads)
	ref   = words[2]
	reads = re.sub('\.', ref, reads)
	reads = re.sub(',', string.lower(ref), reads)

	i = 0
	while i < len(reads):
		i = string.find(reads, '+', i, len(reads))	
		if i < 0:
			break
		match = re.search(r'\d+', reads[i+1:])
                insert_len = int(match.group(0))
		#insert_len = int(reads[i+1])
		reads = reads[:i] + reads[i+insert_len+(match.end()-match.start()) + 1:]

	i = 0
        while i < len(reads):
                i = string.find(reads, '-', i, len(reads))
                if i < 0:
                        break
		match = re.search(r'\d+', reads[i+1:])
                delete_len = int(match.group(0))
                #delete_len = int(reads[i+1])
                reads = reads[:i] + reads[i+delete_len+(match.end()-match.start()) + 1:]

	assert(len(Q) == len(reads))
	qA = ''.join([Q[i] for i, x in enumerate(reads) if x == "A"])
	qT = ''.join([Q[i] for i, x in enumerate(reads) if x == "T"])
	qC = ''.join([Q[i] for i, x in enumerate(reads) if x == "C"])
	qG = ''.join([Q[i] for i, x in enumerate(reads) if x == "G"])

	qa = ''.join([Q[i]  for i, x in enumerate(reads) if x == "a"])
        qt = ''.join([Q[i]  for i, x in enumerate(reads) if x == "t"])
        qc = ''.join([Q[i]  for i, x in enumerate(reads) if x == "c"])
        qg = ''.join([Q[i]  for i, x in enumerate(reads) if x == "g"])		

	fout.write(ref+'\t'+qA+'\t'+qC+'\t'+qG+'\t'+qT+'\t'+qa+'\t'+qc+'\t'+qg+'\t'+qt+'\n')

	
