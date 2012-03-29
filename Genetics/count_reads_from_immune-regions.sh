#!/bin/bash
#$ -cwd

bam=$1

TCRB="7:141998851-142510972"
TCRA="14:22090057-23021075"
IgH="14:106032614-107288051"
MHC="6:29000000-34000000"


x=`samtools view $bam $TCRB | wc -l`
echo -e "TCRB\t$x"


x=`samtools view $bam $TCRA | wc -l`
echo -e "TCRA\t$x"

x=`samtools view $bam $IgH | wc -l`
echo -e "IgH\t$x"

x=`samtools view $bam $MHC | wc -l`
echo -e "MHC\t$x"



