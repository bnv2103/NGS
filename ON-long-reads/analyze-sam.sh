#!/bin/sh
#$ -cwd

prefix=$1

wc -l ${prefix}.sam 
grep Start ${prefix} | cut -f 1 | sort -u | wc -l
grep Start ${prefix} | cut -f 2 | sort -n | uniq -c
grep Start ${prefix} | cut -f 2 | sort -n | uniq -c
grep Start ${prefix} | cut -f 3 | sort | uniq -c
quals=`grep Start ${prefix} | cut -f 5 | sort -n | uniq -c`
echo $quals | tr ' ' '\n' | head -102 | tail -102 | awk 'BEGIN{sum=0;} {if(NR%2==1) sum=sum+$0;} END{print sum;}'
echo $quals | tr ' ' '\n' | head -202 | tail -100 | awk 'BEGIN{sum=0;} {if(NR%2==1) sum=sum+$0;} END{print sum;}'
echo $quals | tr ' ' '\n' | head -302 | tail -100 | awk 'BEGIN{sum=0;} {if(NR%2==1) sum=sum+$0;} END{print sum;}'
echo $quals | tr ' ' '\n' | head -402 | tail -100 | awk 'BEGIN{sum=0;} {if(NR%2==1) sum=sum+$0;} END{print sum;}'
echo $quals | tr ' ' '\n' | head -502 | tail -100 | awk 'BEGIN{sum=0;} {if(NR%2==1) sum=sum+$0;} END{print sum;}' 

