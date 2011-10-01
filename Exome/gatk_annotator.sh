#!/bin/bash
#$ -cwd

heap=4

while getopts v:g:m:b:h opt
  do  
  case "$opt" in
      v) vcf="$OPTARG";;
      g) GLOBAL="$OPTARG";;
      m) MEM="$OPTARG";;
      b) bam="$OPTARG";;
      h) echo $USAGE
	  exit 1;;
  esac
done

if [[ $vcf == "" || $GLOBAL == "" ]]
then
        echo $USAGE
        exit 1
fi

. $GLOBAL

if [[ ! $MEM == "" ]]; then
    heap=$MEM
fi

echo "heap size: ${heap}g"

targets=$vcf".targets.list"

awk '{print $1":"$2"-"$3}' $ExonFile > $targets

TEMP=$vcf"_anno-temp"
mkdir -p $TEMP

JAVA="java -Xmx${heap}g -Djava.io.tmpdir="${TEMP}
GATK="$JAVA -jar "${GATKJAR}

echo $GATK

$GATK \
    -T GenomicAnnotator \
    -R $REF \
    -B:variant,VCF $vcf \
    -L $targets \
    -B:refseq,AnnotatorInputTable $AnnotationTable \
    -o $vcf.genes.annotated \
    -BTI variant \
#    -B:dbsnp,VCF $DBSNPVCF \
 ##   -B:compHapMap,VCF $HapMapV3VCF \
 #   -B:compdbSNP132,VCF $DBSNP132 \
 #   -B:comp1KG,VCF $OneKGenomes
## -m

echo "GenomicAnnotator done."


if [[ $bam != "" ]]; then

    GATK="$JAVA -jar $GATKJAR12"
    $GATK \
	-T VariantAnnotator \
	-R $REF \
	--variant $vcf.genes.annotated  \
	--dbsnp $DBSNPVCF \
	--comp:HapMapV3 $HapMapV3VCF \
	--comp:dbSNP132 $DBSNP132 \
	--comp:1KG $OneKGenomes \
    	-o $vcf.complete.annotated.vcf \
	#	-I $bam \    
    echo "VariantAnnotator done"
fi

rm -rf $TEMP
