import os

p=open('../Regions_1.0_Ion_AmpliSeq_Cancer.bed','r')
p.readline()
d=p.readlines()
p.close()
y=open('output','w')

for item in d:
	nstr=item.split('\t')
	chr=nstr[0]
	start=int(nstr[1])
	end=int(nstr[2])
	q=open('/home/nx2112/dbSNP135/chr_rpts/coor/%s.txt'%chr,'r')
	q.readline()
	q.readline()
	q.readline()
	q.readline()
	q.readline()
	q.readline()
	q.readline()
	t=q.readlines()
	for iter in t:
		new=iter.split('\n')[0]
		if new.isdigit():
			num=int(new)
			if num >=start:
				if num <=end:
					y.write(chr)
					y.write('\t')
					y.write(new)
					y.write('\n')
		else:	
			break
			
	q.close()

y.close()

