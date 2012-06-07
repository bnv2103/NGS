Three types of files are included for each group of data: 

1. summary.csv 
   - contains basic statistical information
	Sample name
	Number of reads
	Number of continually mapped reads
	Median value of FPKM (isoforms)
	Mean value of FPKM (isoforms)
	Number of transcripts with FPKM > 1
	Number of transcripts with FPKM > 0.1
	Number of genes with FPKM > 1
	Number of genes with FPKM > 0.1

2. summary.pdf -
   a.	Histogram of distributions of FPKM in log10 base for all isoforms 
   b.	Plot of top 100 isoforms ranked by normalized FPKM (FPKM *length of isoform)
   c.	(Optional) MA-plots of two replicates: diagnostic plots for assessing the reproducibility of technical or biological replicates.	
   
3. sampleName_isoforms.csv and sampleName_genes.csv 
   - General abundance information of isoforms and genes: gene name, gene id, gene length, FPKM.
 	
4. Differentially expressed gene analysis (optional)
   The following columns are included:
      test_id
      gene_id
      locus
      sample1
      sample2
      status
      value_1: FPKM for sample1
      value_2: FPKM for sample2
      log2.fold: log ratio of fold changes
      test_stat: statistical test score
      p-value
      q_value: adjusted p-value
      significant: identified significant genes
      Notes: genes have low expressed value (FPKM < 1)

5. sampleName_var (optional) 
   variant calling analysis result
   QUAL: Quality score
   DP: Raw Read Depth
   DP4: # high-quality ref-forward bases, ref-reverse bases, alt-forward bases, alt-reverse bases.
   GT: The genotype of this sample
   PL: Phred-scaled likelihoods for AA,AB,BB genotypes where A=ref and B=alt
   GQ: The Genotype Quality

