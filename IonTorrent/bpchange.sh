#!/bin/sh

file=$1

awk '{ if($2 != 4) { for(i=12; i <= NF; i++) {split($i,field,":"); if(field[1]=="MD") {str=field[3];pos=0;num=0;iter=0;pos2=0;num2=0;itt=0;for(it=1;it<=length(str);it++) {sb=substr(str,it,1); if(sb ~ /[0-9]/) {num=10*num+sb;} else if(sb == "^") {pos=pos+num;num=0;it++;while(substr(str,it,1) ~ /[AGCT]/) it++; it--;} else {pos=pos+num+1;num=0; new=sb; iter=0;pos2=0;num2=0;itt=0; while(itt<=length($6)) {itt++; bp=substr($6,itt,1); if(bp ~ /[0-9]/) {num2=10*num2+bp;} else if(bp=="D") {num2=0;} else if(bp ~ /[IS]/) {pos2=pos2+num2;num2=0;} else {iter=iter+num2;pos2=pos2+num2;num2=0;if(iter>=pos) {old=substr($10,pos2-iter+pos,1); print old, "\t", new; break; }}}}}}}}}' $file > bp.txt


grep ^A bp.txt | grep -c C >> AC.txt
grep ^A bp.txt | grep -c G >> AG.txt
grep ^A bp.txt | grep -c T >> AT.txt
grep ^C bp.txt | grep -c G >> CG.txt
grep ^C bp.txt | grep -c T >> CT.txt
grep ^C bp.txt | grep -c A >> CA.txt
grep ^G bp.txt | grep -c T >> GT.txt
grep ^G bp.txt | grep -c A >> GA.txt
grep ^G bp.txt | grep -c C >> GC.txt
grep ^T bp.txt | grep -c A >> TA.txt
grep ^T bp.txt | grep -c C >> TC.txt
grep ^T bp.txt | grep -c G >> TG.txt


