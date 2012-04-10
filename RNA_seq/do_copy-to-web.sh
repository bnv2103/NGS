#!/bin/bash
#$ -cwd

for f in `cat list.recalib_bam.list`; do g=`echo $f | cut -f11 -d '/' | cut -f1 -d '_'`; echo $g;
cp $f /ifs/data/shares/deepsequencing/solid/public_html/CHUNG_12-exomes_2011-09/BAMs/$g
done

