import sys,os,re

#  IUPAC={'a':'A', 'b':'C/G/T','B':'C/G/T','c':'C','d':'A/G/T','D':'A/G/T','g':'G','h':'A/C/T','H':'A/C/T','k':'G/T','K':'G/T','m':'A/C','M':'A/C','n':'A/C/G/T','N':'A/C/G/T','r':'A/G','R':'A/G','s':'C/G','S':'C/G','t':'T','v':'A/C/G','V':'A/C/G','w':'A/T','W':'A/T','y':'C/T','Y':'C/T'}
 
def convertIUPAC(string): 
  IUPAC={'a':'A', 'b':'CGT','B':'CGT','c':'C','d':'AGT','D':'AGT','g':'G','h':'ACT','H':'ACT','k':'GT','K':'GT','m':'AC','M':'AC','n':'ACGT','N':'ACGT','r':'AG','R':'AG','s':'CG','S':'CG','t':'T','v':'ACG','V':'ACG','w':'AT','W':'AT','y':'CT','Y':'CT'}
  D={'A','C','G','T'}
  out=""
  for i,f in enumerate(string.rstrip()):
    if f in D:
      out+="%s" %(f.rstrip())
    elif f in IUPAC:
      out+=IUPAC[f] 
  return out


C={}

for line in open(sys.argv[1]):
  (loc,snp_c)=line.strip().split('\t')[:2]
  C[loc]=snp_c

for lp1 in open(sys.argv[2]):
  (loc_p1,snp_p1)=lp1.strip().split('\t')[:2]
  for lp2 in open(sys.argv[3]):
     (loc_p2,snp_p2)=lp2.strip().split('\t')[:2]
     if loc_p1==loc_p2:
	 if loc_p1 in C:
             print "%s\t%s\t%s" %(convertIUPAC(snp_p1), convertIUPAC(snp_p2), convertIUPAC(C[loc_p1]))
	     break
     else:
         if loc_p1 in C:
             print "%s\t \t%s" %(convertIUPAC(snp_p1), convertIUPAC(C[loc_p1]))
             break
	 else:
           if loc_p2 in C:
	       print " \t%s\t%s" %(convertIUPAC(snp_p2), convertIUPAC(C[loc_p2]))
               break


