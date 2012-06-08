#!/bin/bash
#$ -cwd
#$ -l mem=2G,time=1::

echo "" > $2
grep -w "downstream" $1 | wc -l >> $2
grep -w "upstream" $1 | wc -l >> $2
grep -w "UTR3" $1 | wc -l >> $2
grep -w "UTR5" $1 | wc -l >> $2
grep intergenic $1 | wc -l >> $2
grep intronic $1 | wc -l >> $2
grep ncRNA $1 | wc -l >> $2
grep -w "exonic" $1 | grep -v splicing | wc -l >> $2
grep splicing $1 | wc -l >> $2

