#!/bin/python

# SUMMARY: input a sam file and output data linking quality score to read position and homopolymer proximity


# function homodist
# ----------------
# INPUTS:
#  - minimum homopolymer length
#  - read string

# OUTPUT:
#  - distances homopolymer

def homodist(minhomolen, read):

        homoflags = [0] * len(read)
        # first need to identify homopolymer locations
        posL = 0
        while posL <= len(read) - minhomolen:
                posR = posL
                while posR + 1 < len(read) and read[posR] == read[posR+1]:
                        posR += 1
                homolen = posR - posL + 1
                if homolen >= minhomolen:
                        for homopos in range(posL, posL + homolen):
                                homoflags[homopos] = 1
                posL = posR + 1

	# only looking for homopolymer right boundary
	homoflags[1:] = [int(int(x)-int(y) == -1) for x,y in zip(homoflags[1:], homoflags[0:-1:])]
	homoflags[0] = 0
	homoboundaries = [i-1 for i, e in enumerate(homoflags) if e != 0]
	if homoflags[-1] == 1:
		homoboundaries.append(len(homoflags)-1)
	homodistances = [-1] * len(read)
	if len(homoboundaries) == 0:
		return homodistances 
	for pos in range(len(read)):	
		homodistances[pos] = min([abs(pos - ind) for ind in homoboundaries])
		
        return homodistances


import sys, string

samfilename = sys.argv[1]
samfile = open(samfilename,'r')
lengthoutfile = open(samfilename+".lengths", 'w')
forwardoutfileA = open(samfilename+".forwardA.cor", 'w')
reverseoutfileA = open(samfilename+".reverseA.cor", 'w')
forwardoutfileT = open(samfilename+".forwardT.cor", 'w')
reverseoutfileT = open(samfilename+".reverseT.cor", 'w')
forwardoutfileC = open(samfilename+".forwardC.cor", 'w')
reverseoutfileC = open(samfilename+".reverseC.cor", 'w')
forwardoutfileG = open(samfilename+".forwardG.cor", 'w')
reverseoutfileG = open(samfilename+".reverseG.cor", 'w')

forwardoutfiles = {'A':forwardoutfileA, 'T':forwardoutfileT, 'C':forwardoutfileC, 'G':forwardoutfileG}
reverseoutfiles = {'A':reverseoutfileA, 'T':reverseoutfileT, 'C':reverseoutfileC, 'G':reverseoutfileG}

dic = 'ATCG'

minhomolen = int(sys.argv[2])
assert(minhomolen > 0)
# forwardoutfile.write("# quality\tread_position\thomopolymer_distance_length>="+str(minhomolen)+"\n")
# reverseoutfile.write("# quality\tread_position\thomopolymer_distance_length>="+str(minhomolen)+"\n")

comp = string.maketrans(dic, dic[::-1])

for line in samfile:
        fields = line.split()
	flags = int(fields[1]) # bitwise flag
	if (flags & 4) >> 2:
		continue
	revcomp = (flags & 16) >> 4
        read = fields[9]
	lenread = len(read)
	lengthoutfile.write(str(lenread)+"\n")
        Qstr = fields[10]
	if revcomp:
		read = read[::-1].translate(comp)
		Qstr = Qstr[::-1]
        homodistances = homodist(minhomolen, read)
        for i in range(lenread):
		if revcomp:
                	reverseoutfiles[read[i]].write(str(ord(Qstr[i]))+"\t"+str(i)+"\t"+str(homodistances[i])+"\t"+str(lenread)+"\n")
		else:
			forwardoutfiles[read[i]].write(str(ord(Qstr[i]))+"\t"+str(i)+"\t"+str(homodistances[i])+"\t"+str(lenread)+"\n")

lengthoutfile.close()
for letter in dic:
	forwardoutfiles[letter].close()
	reverseoutfiles[letter].close()
