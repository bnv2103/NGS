#!/bin/sh
# Code to plot error rate against read length
# Input is file containing positions of all errors to be plotted.
# > sh plot.sh <file> <SAM file>

file=$1
sam_file=$2

sort -n $file | uniq -c > errors.txt
awk '{print length($10)}' $sam_file> | sort -rn | uniq -c > readlen.txt
# Problem arises when #lines in errors.txt and readlen.txt don't match.
# The join command below will fail.
# Need to do a little dirty work here to append "missing" lines, or lines
# pertaining to positions that contain no errors. Sample line to be appended:
#       200 0
# 200 being position, 0 being # errors. Needs to be automated ideally.
awk 'BEGIN{sum=0;} { sub(/^[ \t]+/, ""); sum=sum+$1; print sum, "\t", $2; }' readlen.txt > readlens.txt
awk '{ sub(/^[ \t]+/, ""); print }' errors.txt > error.txt
sort -k2 -nr error.txt | tr ' ' '\t' > errors.txt
sort -k2 -n errors.txt | awk '{print $2, "\t", $1}' > err.txt
sort -k2 -n readlens.txt | awk '{print $2, "\t", $1}' > rd.txt
join -1 1 err.txt rd.txt | tr ' ' '\t' | awk '{ print $1, "\t", $2/$3}' > plot.txt
gnuplot
# set term png
# set outout 'output_file'
# plot 'plot.txt'

