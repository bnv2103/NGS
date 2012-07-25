#!/bin/bash
#$ -cwd
TEMP=$1"_temp1"
if [ ! -e $TEMP ]
then
mkdir $TEMP
fi

JAVA="java -Xmx2G -Djava.io.tmpdir="${TEMP}
PJAR="/ifs/data/c2b2/ngs_lab/ngs/usr/picard-tools-1.65/CreateSequenceDictionary.jar"

outf=`echo $1 | sed 's/.fasta/.dict/' `
$JAVA -jar $PJAR R="$1" O="$outf"
