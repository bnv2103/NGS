Three types of files are included for each group of data: 

1. summary.csv 
   - contains basic statistical information
	Sample name
	Number of raw reads
	Number of mapped reads
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
   - General information of isoforms and genes: gene name, gene id, gene length, FPKM.
 	

