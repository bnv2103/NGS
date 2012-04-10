#!/bin/bash
#$ -cwd

f1=$1
f2=$2

# get unique reads
#samtools view $f1 | cut -f1 | sort -u -S 2G > readsName1.txt
#samtools view $f2 | cut -f1 | sort -u -S 2G > readsName2.txt
samtools view -X $f1 | cut -f1,2 | grep 'pP' | cut -f1 |sort -u -S 18G > readsName1PE.txt
samtools view -X $f2 | cut -f1,2 | grep 'pP' | cut -f1 |sort -u -S 18G > readsName2PE.txt

# get overlaped reads
comm -12 readsName1PE.txt readsName2PE.txt > commReadsPE.txt
wc -l readsName1PE.txt > numOfSample1
wc -l readsName2PE.txt > numofSample2
wc -l commReadsPE.txt > numOfComm
# get reads name and MAPQ value
#samtools view $f1 | cut -f1,5 > reads1.txt
#samtools view $f2 | cut -f1,5 > reads2.txt

# sort reads based on their MAPQ 
#for i in `cat commReads.txt`; do grep -w $i reads1.txt |cut -f2 |sort |tail -1> a1;
#    grep  -w $i reads2.txt |cut -f2 |sort |tail -1> a2;
#   if [ `cat a1` -gt `cat a2` ]; then
#       grep -w $i reads2.txt | cut -f1 |sort -u  >> results4Ms.txt
#   elif [ `cat a2` -gt `cat a1` ]; then
#       grep -w $i reads2.txt | cut -f1 |sort -u  >> results4Hm.txt
#   else
#       grep -w $i reads2.txt | cut -f1 |sort -u  >> results4Both.txt
#   fi
#done

#wc -l results4Ms.txt > MsCount
#wc -l results4Hm.txt > HmCount
#wc -l results4Both.txt > both
