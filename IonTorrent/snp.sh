#!/bin/sh

# file is SAM file>
file=$1

# Algorithm for the so easily decipherable awk code below
#
# starting from 12th field until the end
#   split on : and check for MD in first part
#     if(MD) obtain third part
#       while(parse)
#         if(number) num=10*num+number
#         if(^) skip all characters ahead and pos=pos+num
# 	  while(char) echo pos
#	end parse
#     end MD
#   end split
# end

awk '{ if($2 != 4) { for(i=12; i <= NF; i++) {split($i,field,":"); if(field[1]=="MD") {str=field[3];pos=0;num=0;for(it=1;it<=length(str);it++) {sb=substr(str,it,1); if(sb ~ /[0-9]/) {num=10*num+sb;} else if(sb == "^") {pos=pos+num;num=0;it++;while(substr(str,it,1) ~ /[AGCT]/) it++; it--;} else {pos=pos+num;num=0;if($2==16) {rev=length($10);print rev-pos+1;} else print pos; }}}}}}' $file

