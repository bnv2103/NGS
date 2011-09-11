###set global config
export REFTYPE="build"  ## build: no "chr" in chromosome names, aka >1 >2 etc;  hg: >chr1 >chr2 etc.
export REF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/human_g1k_v37.fasta"
#export ExonFile="/ifs/data/c2b2/af_lab/saec/Sequencing/resources/exomes/agilent_1.1_refseq_plus_3_boosters_hg19/agilent_1.1_refseq_plus_3_boosters_b37"
export ExonFile="/ifs/data/c2b2/ngs_lab/ngs/resources/Agilent/SureSelect_All_Exon_50mb_with_annotation.hg19.bed.mod"
export DBSNP="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/dbsnp_130_b37.rod"
export DBSNPVCF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/dbsnp_132.b37.excluding_sites_after_129.vcf" 
export HapMapV3VCF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/hapmap_3.3.b37.vcf"
export INDELVCF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/1000G_indels_for_realignment.b37.vcf"
export SAMTOOLS="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/samtools"
# export TEMP="/ifs/scratch/c2b2/af_lab/saec/temp/"  # temp dir for Java
export GATKJAR="/ifs/data/c2b2/ngs_lab/ngs/usr/GATK/Sting/dist/GenomeAnalysisTK.jar" 
export PICARD="/ifs/data/c2b2/ngs_lab/ngs/usr/picard-tools"
export BWA="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/bwa"  
export BPATH="/ifs/home/c2b2/ngs_lab/ngs/code/NGS/Exome"
export STING="/ifs/data/c2b2/ngs_lab/ngs/usr/GATK"
export AnnotationTable="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/refGene-big-table-b37.txt"

