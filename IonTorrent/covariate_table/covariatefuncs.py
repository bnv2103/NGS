#!/bin/python

from re import compile, finditer, findall

def homodist(read):

    matchobjs = list(compile('(A{3,}|T{3,}|G{3,}|C{3,})').finditer(read))
    homoboundaries = [x.end() - 1 for x in matchobjs]
    if not homoboundaries:
        return [0]*len(read) 

    out = [-1]*len(read)

    matches = [x.group() for x in matchobjs]

    for pos in range(len(read)):
        Ds = [abs(pos - b) for b in homoboundaries]
        minD = min(Ds)
        ind = Ds.index(minD) 
        #if homoboundaries[ind] < pos:
        #    homodistance = minD - 1
        #else:
        homodistance = minD
        if homodistance > 1 or (homodistance > 0 and homoboundaries[ind] > pos):
            out[pos] = 0
        else:
            out[pos] = [len(matches[ind]), homodistance]

    return out 

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

