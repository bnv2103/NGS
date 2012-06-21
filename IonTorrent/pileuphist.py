#!/bin/python

import sys

fin = open(sys.argv[1], 'r')
fout = open("pileuphist.out", 'w')

fout.write("position\tdepth\tmismatch\tref\tallele1\tallele2\tallele3\n")

for line in fin:
	words = line.split()
	reads = words[4] 
	depth = words[3] 
	pos   = words[1]
	ref   = words[2]
	ref_i = "ATCG".find(ref)
	mismatch = 0
	readi = 0
	counts = [0, 0, 0, 0]
	while readi < len(reads):
		read = reads[readi]
		if read in ".,":
			readi = readi + 1
			###ref_i = "ATCG".find(ref)
			counts[ref_i] = counts[ref_i] + 1
		elif read in "ACTGNactgn":
			readi = readi + 1
			mismatch = mismatch + 1 
			if read in "Aa":
				counts[0] = counts[0] + 1
			elif read in "Tt":
				counts[1] = counts[1] + 1
			elif read in "Cc":
                                counts[2] = counts[2] + 1 
			elif read in "Gg":
                                counts[3] = counts[3] + 1 
		elif read == '+':
			insert_length = int(reads[readi+1])
			readi = readi + insert_length + 2
			mismatch = mismatch + 1
		elif read == '-':
                        delete_length = int(reads[readi+1])
                        readi = readi + delete_length + 2
                        mismatch = mismatch + delete_length 
		elif read == '^':
			readi = readi + 2
		elif read == '$':
			readi = readi + 1
		elif read == '*':
			readi = readi + 1
		else:
			print "Error: unrecognized read character \""+read+'\"'+" in read string \""+reads+'\"'
			sys.exit()
	
	ref_ct = counts.pop(ref_i)
	allele_cts = sorted(counts, reverse=True)
	fout.write(str(pos)+'\t'+str(depth)+'\t'+str(mismatch)+'\t'+str(ref_ct)+'\t'+str(allele_cts[0])+'\t'+str(allele_cts[1])+'\t'+str(allele_cts[2])+'\n')

	
