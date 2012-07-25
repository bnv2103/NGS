job.sh is the master script. Any parameters need to be specified in this script.
Presently indel calling is failing with Java Heap out of space error. I've allocated 1G of memory for the Heap. Allocating more than 1G of memory for the Heap is failing. Tried testing using qsub but the job is never triggered. Hence I am not able to test beyond indel calling. 

##################
Date: July 24, 2012
Samreen Zafer 
ToolUpdate: Gatk has been updated to version 1.16.13
ToolUpdate: Annovar has been update to Version: $LastChangedDate: 2012-05-25 01:57:38 -0700 (Fri, 25 May 2012) and convert2annovar.pl updated on july 18,2012

-All Exome Seq and Annovar scripts have been update to use latest GATK parameters.
-Auto pipeline now calls Annovar for annotation and then does filtering.
-gatk_combine_variants will now merge all the input VCFs from diff samples and then perform Annovar->annotation followed by filtering.
-Gatks VQSR has been implemented. I plan to replace the hard filtering on SNV (gatk_filter.sh) with gatk_VQSR_snp_indel_WES.sh , which is then followed by vcf_filter-indel.sh

OlderScripts now reside in directory ~NGS/Exome_Gatk1.0/
##################
