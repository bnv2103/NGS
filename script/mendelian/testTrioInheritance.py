import sys,os,math

# Set track values
trialleles = 0;
parentallele = 0;
nexamined = 0;
ncorrect = 0;
nerrors = 0;

for line in open(sys.argv[1]):
  (parent1splt, parent2splt, childsplt)=line.rstrip().split('\t')[:3]
  if len(childsplt)>2:
    trialleles += 1
  else:
    if (len(parent1splt)>2 and len(parent2splt)>2):
      if len(childsplt)>=2:
        parentallele += 1
      else:
        if (len(parent1splt)>2 or len(parent2splt)<=2) :
          if parent2splt.find(childsplt):
            ncorrect += 1
          else:
            parentallele += 1
        if (len(parent1splt)<=2 or len(parent2splt)>2) :
          if parent1splt.find(childsplt) :
            ncorrect += 1
          else:
            parentallele += 1
        if (len(parent1splt)>2 or len(parent2splt)>2) :
            parentallele += 1
    else:
      if len(childsplt)==2:
        if (parent1splt.find(childsplt) and parent2splt.find(childsplt)) :
          if parent1splt.find(childsplt):
            if parent2splt.find(childsplt):
              ncorrect += 1
          if (ncorrect==0 or parent2splt.find(childsplt)) :
            if parent1splt.find(childsplt):
              ncorrect += 1
        if ncorrect == 0:
          nerrors += 1
      else:
        if (parent1splt.find(childsplt) and parent2splt.find(childsplt)) :
          ncorrect += 1
        else:
          nerrors += 1
  nexamined+=1
print " Number of observed SNPs in child: %s \n Number of SNPs attributable to inheritance: %s \n Number of SNPs that represent Mendelian errors: %s \n Number of SNPs with unclear parent calls (excluded): %s \n Number of tri-allelic SNPs in child (excluded): %s " %(nexamined,ncorrect,nerrors,parentallele,trialleles)
falsePositive=float(nerrors)/(float(nerrors)+float(ncorrect))
print " False Positive Error Rate: %f" %(falsePositive)

