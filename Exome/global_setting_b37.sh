###set global config
export REFTYPE="build"  ## build: no "chr" in chromosome names, aka >1 >2 etc;  hg: >chr1 >chr2 etc.
export REF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/human_g1k_v37.fasta"
#export ExonFile="/ifs/data/c2b2/af_lab/saec/Sequencing/resources/exomes/agilent_1.1_refseq_plus_3_boosters_hg19/agilent_1.1_refseq_plus_3_boosters_b37"
export ExonFile="/ifs/data/c2b2/ngs_lab/ngs/resources/Agilent/SureSelect_All_Exon_50mb_with_annotation.hg19.bed.mod"
export DBSNP="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/dbsnp_130_b37.rod"
export DBSNPVCF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/dbsnp_132.b37.excluding_sites_after_129.vcf" 
export DBSNP132="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/dbsnp_132.b37.vcf"
export DBSNP135="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/dbsnp_135_b37.vcf"
export HapMapV3VCF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/hapmap_3.3.b37.vcf"
export INDELVCF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/1000G_indels_for_realignment.b37.vcf"
export OneKGenomes="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/ALL.wgs.phase1.projectConsensus.snps.sites.vcf.PASS.vcf"
export EVSVCF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/EVS.vcf"
export SAMTOOLS="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/samtools"
# export TEMP="/ifs/scratch/c2b2/af_lab/saec/temp/"  # temp dir for Java
export GATKJAR_old="/ifs/data/c2b2/ngs_lab/ngs/usr/GATK/Sting/dist/GenomeAnalysisTK.jar" 
export GATKJAR12="/ifs/data/c2b2/ngs_lab/ngs/usr/GATK-1.2-1/GenomeAnalysisTK.jar"
export GATKJAR="/ifs/scratch/c2b2/ngs_lab/sz2317/softwares/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar"
export PICARD="/ifs/data/c2b2/ngs_lab/ngs/usr/picard-tools"
export BWA="/ifs/data/c2b2/ngs_lab/ngs/usr/bin/bwa"  
export BPATH="/ifs/home/c2b2/ngs_lab/ngs/code/NGS/Exome"
export STING="/ifs/data/c2b2/ngs_lab/ngs/usr/GATK"
export AnnotationTable="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/refGene-big-table-b37.txt"
export RUBY18="/ifs/data/c2b2/ngs_lab/ngs/usr/local/bin/ruby"
export UTILS="/ifs/data/c2b2/ngs_lab/ngs/code/NGS/utils/"
export ExonFile=/ifs/data/c2b2/ngs_lab/ngs/resources/Agilent/SureSelect_All_Exon_V4_hg19.bed.sorted.bed
export ANNOVAR="/ifs/scratch/c2b2/ngs_lab/sz2317/softwares/annovar/"
export OMNI_Sites="/ifs/scratch/c2b2/ngs_lab/sz2317/softwares/GATK_resourceBundle/1000G_omni2.5.b37.sites.vcf"
export HamMap_Sites="/ifs/scratch/c2b2/ngs_lab/sz2317/softwares/GATK_resourceBundle/hapmap_3.3.b37.sites.vcf"
export Mills1KG_Sites="/ifs/scratch/c2b2/ngs_lab/sz2317/softwares/GATK_resourceBundle/Mills_and_1000G_gold_standard.indels.b37.sites.vcf"

