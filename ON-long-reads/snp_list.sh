#!/bin/sh
#$ -cwd

snp=$1

#grep -v ^@ test_*.sam | sort -k4,4n | awk -v snp=$snp 'BEGIN{ s[1]=9888016;s[2]=9893844;s[3]=9894822;s[4]=9917426;s[5] = 9919921 }

grep -v ^@ test_10.sam | sort -k4,4n | awk -v snp=$snp 'BEGIN{ s[1]=9883185;s[2]=9887958;s[3]=9888048;s[4]=9893648;s[5]=9903647;s[6]=9906538;s[7]=9917998;s[8]=9918151;s[9]=9919921;s[10]=9927942; }
	{ start = $4-1; cigar = $6;
	for(i=snp;i<=snp;i++) {
		actual = 0; num = 0; pos = 0;
		for(it=1;it<=length(cigar);it++) {
			sb = substr(cigar,it,1);
			if(sb ~ /[0-9]/) { num=num*10+sb;}
			else if(sb ~ /[MS]/) { actual+=num;pos+=num;num=0;}
			else if(sb == "D") {actual+=num; num=0; }
			else if(sb == "I") {pos+=num; num=0; }
			else print "Encountered", sb;

			if(actual+start >= s[i] && start <= s[i]) {
				##print NR,actual,pos,s[i],start;
				pos = pos-actual+s[i]-start;
				print substr($10,pos,1), $1;
				break;
			}
		}
	}
}' 

