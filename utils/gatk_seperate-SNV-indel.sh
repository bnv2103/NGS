#!/bin/sh
#$ -S /bin/sh
#$ -cwd
uname -a

# Filter SNPs in indel regions
# Load all samples from chromosome in a single $INP vcf file


HEAP=2
INP=""
REF=""

USAGE="Usage: $0 -I <Input VCF file> -g "

while getopts I:t:m:g:h:A: o
do      case "$o" in
        I)      INP="$OPTARG";;
        m)      MEM="$OPTARG";;
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

# JOB_ID is the qsub job ID
TEMP=$INP"_sel"

if [ ! -d $TEMP ]; 
    then    
	mkdir -p $TEMP
fi

JAVA="java -Xmx${HEAP}g -Djava.io.tmpdir="${TEMP}
GATK="$JAVA -jar "${GATKJAR}

$GATK \
   -R $REF \
   -T SelectVariants \
   --variant $INP \
   -o $INP.snv \
   -selectType SNP 

echo "SNV select complete"

$GATK \
   -R $REF \
   -T SelectVariants \
   --variant $INP \
   -o $INP.indel \
   -selectType INDEL 

echo "Indel select complete"

