#!/bin/bash
#$ -cwd

HEAP=4
threads=1

USAGE="Usage: $0 -v vcf -g global_settings -b bam.list [ -o output] [ -m heap_size] [-n threads]"

while getopts v:b:o:g:n:m:h opt
  do  
  case "$opt" in
      v) input="$OPTARG";;
      b) bam="$OPTARG";;
      g) settings="$OPTARG";;
      o) output="$OPTARG";;
      n) threads="$OPTARG";;
      m) HEAP="$OPTARG";;
      h)      echo $USAGE
	  exit 1;;
  esac
done

if [[ $input == "" || $settings == "" || $bam == "" ]]
    then
    echo $USAGE
    exit 1
fi

. $settings 

temp=$input"_annot_temp"
mkdir -p $temp

if [[ $output == "" ]] 
    then output=$input".moreinfo.vcf"
fi


JAVA="java -Xmx${HEAP}g -Djava.io.tmpdir="${temp}
GATK="$JAVA -jar "${GATKJAR12}

uname -a 

$GATK \
    -T VariantAnnotator \
    -R $REF \
    -o $output \
    --variant $input \
    --dbsnp $DBSNPVCF \
    --comp:HapMapV3 $HapMapV3VCF \
    --comp:dbSNP132 $DBSNP132 \
    --comp:1KG $OneKGenomes 
#    -all \
#    -I $bam \
#    -B:comp1KG_CEU,VCF 1000GenomesCalls.CEU.vcf \

#    -A AlleleBalance \
#    -A DepthOfCoverage \
#    -A BaseQualityRankSumTest \
#    -A HomopolymerRun \
#    -A MappingQualityRankSumTest \
#    -A MappingQualityZero \
#    -A QualByDepth \
#    -A RMSMappingQuality \
#    -A SpanningDeletions \
#    -A HaplotypeScore  \


rm -rf $temp
