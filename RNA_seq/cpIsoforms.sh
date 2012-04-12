#!/bin/bash
#$ -cwd

infile=$1
outdir=$2

# copy isoforms and genes information generated from cufflink to current working directory

isoforms=$infile"/accepted_hits.bam_cufflinks_ref/isoforms.fpkm_tracking"
genes=$infile"/accepted_hits.bam_cufflinks_ref/genes.fpkm_tracking"
SNPs=$infile"/SNPs/var_summary.txt"

infileB=` echo $infile | cut -f1 -d '.'`
cp $isoforms $outdir"/"$infileB"_isoforms"
cp $genes $outdir"/"$infileB"_genes"

infileC=` echo $infile | cut -f6 -d '_'`
cp $SNPs $outdir"/"$infileC"_var"

# echo $infileB
 


