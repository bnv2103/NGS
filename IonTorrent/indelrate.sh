#!/bin/sh
# read size vs overall error rate of a read

file=$1

awk '{if($2!=4){err=0;for(i=1;i<length($6);i++) {if(substr($6,i,1)=="I"||substr($6,i,1)=="D") err++;} print length($10), "\t", err;}}' $file | sort -k1,1n | awk 'BEGIN{readlen=0;ct=0;sum=0;} {if($1!=readlen) {if(readlen!=0){print readlen, "\t", sum/(ct);}readlen=$1;sum=0;ct=0;} else {ct++;sum=sum+$2;}} END{print readlen, "\t", sum/(ct);}' > indel_rate.txt
awk '{if($2!=4){err=0;for(i=12;i<=NF;i++){split($i,field,":"); if(field[1]=="MD") {str=field[3];for(it=1;it<=length(str);it++) {sb=substr(str,it,1); if(sb ~ /[0-9]/) continue; else if(sb == "^") {it++;while(substr(str,it,1) ~ /[AGCT]/) it++; it--;} else err++;}}} print length($10), "\t", err;}}' $file | sort -k1,1n | awk 'BEGIN{readlen=0;ct=0;sum=0;} {if($1!=readlen) {if(readlen!=0){print readlen, "\t", sum/(ct);}readlen=$1;sum=0;ct=0;} else {ct++;sum=sum+$2;}} END{print readlen, "\t", sum/(ct);}' > snp_rate.txt

