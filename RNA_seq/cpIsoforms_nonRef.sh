#!/bin/bash
#$ -cwd

infile=$1
outdir=$2

# copy isoforms and genes information generated from cufflink to current working directory

isoforms=$infile"/accepted_hits.bam_cufflinks_ref-guide/isoforms.fpkm_tracking"
genes=$infile"/accepted_hits.bam_cufflinks_ref-guide/genes.fpkm_tracking"
infileB=` echo $infile | cut -f1 -d '.'`
cp $isoforms $outdir"/"$infileB"_isoforms_nonRef"
cp $genes $outdir"/"$infileB"_genes_nonRef"
# echo $infileB
 


