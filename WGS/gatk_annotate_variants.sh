#!/bin/bash
#$ -cwd
#$ -l mem=8G,time=1::

uname -a

REF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/human_g1k_v37.fasta"
GATK="/ifs/scratch/c2b2/ngs_lab/sz2317/softwares/GenomeAnalysisTK-1.5-21-g979a84a/GenomeAnalysisTK.jar"
# VCF="/ifs/scratch/c2b2/ngs_lab/ngs/Projects/WG-seq/120309_ANDREA_GABRIELLE_6_HUMAN_GENOME_30X_PE_HISEQ/SAMPLES/AC2/VarCalling_samtool/1.var.flt.vcf"
 VCF=$1

DBSNP135="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/dbsnp_135_b37.vcf"
HapMapV3VCF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/hapmap_3.3.b37.vcf"
OneKGenomes="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/ALL.wgs.phase1.projectConsensus.snps.sites.vcf.PASS.vcf"
EVSVCF="/ifs/data/c2b2/ngs_lab/ngs/resources/bwa_samtools_gatk_DB/EVS.vcf"

TEMP=$VCF"_temp1"
if [ ! -e $TEMP ]
then
mkdir $TEMP
fi

# List of annotations from /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/Exome/Exome_newGATK/Variant_Annotation_List.txt.new 

java -Xmx3g -Djava.io.tmpdir=$TEMP -jar $GATK \
	-R $REF \
	-I $2 \
	-T VariantAnnotator  \
	--variant $VCF \
	-o $VCF".annotate_varnt_all.vcf" \
	-L $VCF  \
	-XA "SnpEff" \
	-XA "TransmissionDisequilibriumTest" \
	-XA "ChromosomeCounts" \
	-XA "HardyWeinberg" \
	-XA "NBaseCount" \
	-XA "BaseCounts" \
	-XA "MVLikelihoodRatio" \
	-XA "RodRequiringAnnotation" \
	-XA "TechnologyComposition" \
	-XA "SampleList" \
	-all \
        -A AlleleBalance  -A DepthOfCoverage  -A BaseQualityRankSumTest  -A HomopolymerRun -A MappingQualityRankSumTest -A MappingQualityZero -A QualByDepth  -A RMSMappingQuality  -A SpanningDeletions  -A HaplotypeScore


#Do ANNOVAR genomic annotation
echo " qsub -o $1.log.annovar.o -e $1.log.annovar.e /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/do_annovar_all.sh $VCF.annotate_varnt_all.vcf "
#  qsub -o $1.log.annovar.o -e $1.log.annovar.e /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/do_annovar_all.sh $VCF.annotate_varnt_all.vcf 


