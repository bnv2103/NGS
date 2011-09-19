#!/bin/bash
#$ -cwd

vcf=$1

head -600 $vcf | egrep "^#" >  $vcf.header
egrep -v "^#" $vcf  | grep -w PASS > $vcf.pass
cat $vcf.header $vcf.pass > $vcf.release
