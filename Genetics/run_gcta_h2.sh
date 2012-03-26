#!/bin/bash
#$ -cwd

prefix=$1
out=$2
prev=$3

if [[ $prev == "" ]]; then
    prev=0.0001
fi

gcta --bfile $prefix --autosome --maf 0.01 --make-grm --out $out.grm

awk '{print $1,$2,$6}' $prefix.fam > $out.pheno.txt
gcta64 --grm $out.grm --pheno $out.pheno.txt --reml --prevalence $prev --out $out.h2


