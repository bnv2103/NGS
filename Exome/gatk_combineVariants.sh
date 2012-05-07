#!/bin/sh
#$ -S /bin/sh
#$ -cwd
uname -a

# iCombine Variants from different samples (each in one vcf file)

HEAP=2
INPLIST=""
REF=""

USAGE="Usage: $0 -L <File with list of Input VCF files> -g "

while getopts L:O:t:m:g:h:A: o
do      case "$o" in
        L)      INPLIST="$OPTARG";;
        m)      MEM="$OPTARG";;
	O)	OUTVCF="$OPTARG";;
        g)      GLOBAL="$OPTARG";;
	A)	AUTO="$OPTARG";;
        h)      echo $USAGE
                exit 1;;
        esac
done

if [[ $INPLIST == "" || $GLOBAL == "" ]]
    then
    echo $USAGE
    exit 1
fi

. $GLOBAL

if [[ $OUTVCF == "" ]]; then OUTVCF="combined"; fi
if [ ! $MEM == "" ]
    then
    HEAP=$MEM
fi

dname=`dirname $INPLIST`
d1=`dirname $dname`
sname=`basename $d1`
if [ $AUTO == "" ]; then
	job_ext="$sname"
else
	job_ext="$sname.AUTO"
fi
if [ $JOB_ID == "" ]; then
    JOB_ID="VarFilter$job_ext"
fi

cd $dname
# JOB_ID is the qsub job ID
TEMP=$INPLIST"_"$JOB_ID"/"

if [ ! -d $TEMP ]; 
    then    
	mkdir -p $TEMP
fi

if [ -e $dname"/*.idx" ];
then
	rm $dname"/*.idx"
fi

JAVA="java -Xmx${HEAP}g -Djava.io.tmpdir="${TEMP}
GATK="$JAVA -jar "${GATKJAR15}

#  $GATK -T CombineVariants -R $REF -I $INPLIST  -o $OUTVCF.vcf -genotypeMergeOptions UNSORTED -filteredRecordsMergeType KEEP_UNCONDITIONAL
cmd=" $GATK -T CombineVariants -R $REF "
 for i in `cat $INPLIST`;do 
	cmd="$cmd --variant $i "
 done
 cmd="$cmd -o $OUTVCF.vcf -genotypeMergeOptions UNSORTED -filteredRecordsMergeType KEEP_UNCONDITIONAL " 

echo $cmd
$cmd

