#!/bin/bash
#$ -cwd

SAMTOOLS=`which samtools`

if [[  $SAMTOOLS == ""  ]]; then
    SAMTOOLS="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/samtools"
fi


${SAMTOOLS} index $1 
echo "$1 indexed"
