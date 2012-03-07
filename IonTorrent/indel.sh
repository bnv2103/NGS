#!/bin/sh

# Code to extract position of all indels from a given SAM file for identifying
# error patterns of a sequencing machine.
#
# To run, > sh indel.sh <SAM file>  >  <OP file>
#
# This will return the set of positions where indels exist, which can then be
# redirected to another file for plotting.
# This is not efficient code. It will run much faster using awk (as
# implemented for snps already).

file=$1

awk '{if($2!=4){cigar=$6;num=0;pos=0; for(i=1;i<=length(cigar);i++){bp=substr(cigar,i,1);if(bp ~ /[0-9]/) {num=10*num+bp;} else if(bp ~ /[MS]/) {pos=pos+num;num=0;} else {if(bp=="I"){pos=pos+num;}num=0; if($2==16) {print length($10)-pos+1;} else print pos; }}}}' $file

