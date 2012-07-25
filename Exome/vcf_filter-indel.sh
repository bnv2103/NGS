#!/bin/bash
#$ -cwd
uname -a

# Filter indels
# Load all samples from chromosome in a single $INP vcf file


HEAP=4
INP=""
REF=""

USAGE="Usage: $0 -I <Input VCF file> -g "

while getopts I:t:m:g:h:A:s: o
do      case "$o" in
        I)      INP="$OPTARG";;
        m)      MEM="$OPTARG";;
	s)	SNV="$OPTARG";;
        g)      GLOBAL="$OPTARG";;
	A)	AUTO="$OPTARG";;
        h)      echo $USAGE
                exit 1;;
        esac
done

if [[ $INP == "" || $GLOBAL == "" ]]
    then
    echo $USAGE
    exit 1
fi

. $GLOBAL

if [ ! $MEM == "" ]
    then
    HEAP=$MEM
fi

if [ $JOB_ID == "" ]; then
    JOB_ID="IndelFilter"
fi

# JOB_ID is the qsub job ID
TEMP=$INP"_"$JOB_ID"/"

if [ ! -d $TEMP ]; 
    then    
    mkdir -p $TEMP
fi

JAVA="java -Xmx${HEAP}g -Djava.io.tmpdir="${TEMP}
GATK="$JAVA -jar "${GATKJAR}

$GATK \
    -T VariantFiltration \
    -R $REF \
    -o $INP.filtered.vcf \
    --variant $INP \
    --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)"  \
    --filterName "HARD_TO_VALIDATE" \
    --filterExpression "SB >= -1.0" \
    --filterName "StrandBiasFilter" \
    --filterExpression "QUAL < 10" \
    --filterName "QualFilter"

sh ${BPATH}/do_release.sh $INP.filtered.vcf

if [[ ! $SNV == "" ]];then  ## Merge the filtered snv and indel variants together.
	CMD="sh $UTILS/gatk_merge_SNP_Indel.sh -s $SNV -i $INP.filtered.vcf -o $INP.allVariants.filtered.vcf "
	echo $CMD
	$CMD
fi

