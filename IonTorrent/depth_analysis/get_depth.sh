#!/bin/bash
# -cwd

for s in `ls ../sample_files/*.bam`
do
    sample=$s.mpileups
    for g in `ls ../sample_files/$sample`
    do
        for f in `ls ../sample_files/$sample/$g`
        do
             echo -e "`basename $f .mpileup`\t`awk -F'\t' '{ sum += $4 } END { print sum / NR }' ../sample_files/$sample/$g/$f`" >> `basename $s .bam`.depth 
        done
    done
done
