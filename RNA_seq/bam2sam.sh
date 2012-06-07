#!/bin/bash
#$ -cwd

# gtf_mouse="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf"
#  gtf_human="/ifs/data/c2b2/ngs_lab/ngs/resources/bowtie_DB/Homo_sapiens/NCBI/build37.2/Annotation/Genes/genes.gtf"



outDir="DESEQ/"
mkdir $outDir

for f in *cufflinks

do
 samFile=$f".sam"
 txtFile=$f".txt"
 samtools view -h -o $outDir$samFile $f"/accepted_hits.bam"
 
#  htseq-count $outDir$samFile $gtf_human --stranded=no > $outDir$txtFile
done
