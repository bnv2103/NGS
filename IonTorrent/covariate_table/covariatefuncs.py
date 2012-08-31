#!/bin/python

from re import compile, finditer, findall

# this function takes a sequence of nucleotides as input and outputs a list of the homopolymer info at each position on the sequence
# this info (each element of the list) is either 0 (signifying no homopolymer to worry about) or it is a list of two integers; the 
# length of the homopolymer that is nearby, and how near its 3' end is from the current position. We only consider distance of zero or
# one base downstream from the 3' end as relevant, and only homopolymer of length 3 or greater (other wise the position gets the null 0 entry)
def homodist(read):

    # find the homopolymers of at least length 3
    matchobjs = list(compile('(A{3,}|T{3,}|G{3,}|C{3,})').finditer(read))
    # get the 3' boundaries for each homopolymer
    homoboundaries = [x.end() - 1 for x in matchobjs]
    # if there are no homopolyers of length 3, return a list of all zeros
    if not homoboundaries:
        return [0]*len(read) 

    # initialize the output list to all -1
    out = [-1]*len(read)

    # list of the explicit homopolymer strings
    matches = [x.group() for x in matchobjs]

    # loop through the nucleotides in the read and compute the corresponding homopolymer output element
    # there should be no more -1 in the out list
    for pos in range(len(read)):
        # distances of all homopolymer 3' boundaries from this position
        Ds = [abs(pos - b) for b in homoboundaries]
        # the closest one, and its index
        minD = min(Ds)
        ind = Ds.index(minD) 

        homodistance = minD
        # if it's further not on the 3' boundary, or one downstream, forget it
        if homodistance > 1 or (homodistance > 0 and homoboundaries[ind] > pos):
            out[pos] = 0
        else:
            out[pos] = [len(matches[ind]), homodistance]

    return out 

# this function takes a cigar string and a read and a quality string, and outputs the read and quality in alignment form
# NOTE: this function is currently not used, which is good because it is buggy (compare with Matlab's cigar2align results)
def cigar2align(cigar, read, quality):
    read = list(read)
    quality = list(quality)
    cigar = findall('\d+|[M I D N S H P]', cigar)
    idx = 0
    for i in range(0, len(cigar), 2):
        length = int(cigar[i])
        operation = cigar[i+1]
        #if operation == 'H':
        if operation == 'S':
            read = read[:idx] + read[idx + length:]
            quality = quality[:idx] + quality[idx + length:]
        elif operation == 'I':
            read = read[:idx] + read[idx + length:]
            quality = quality[:idx] + quality[idx + length:]
        elif operation == 'D':
            read = read[:idx] + list('-'*length) + read[idx:]
            quality = quality[:idx] + list('-'*length) + quality[idx:]
            idx = idx + length
        elif operation == 'M':
            idx = idx + length
        elif operation not in ['H', 'M']:
            print 'unrecognized cigar operation'
            sys.exit()
    return ''.join(read), ''.join(quality)

