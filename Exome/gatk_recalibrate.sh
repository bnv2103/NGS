#!/bin/sh
#$ -cwd
# Findmem

# Recalibrate base quality score
# INPUT:	Run with a single read-group/lane in $INP.list
# OUTPUT:	$INP.$CHR.recalibrated.bam recalibrated bam file

HEAP=4000

# CHRFILE="${BPATH}/chrdb"

INP=""
CHR=""
REF=""
OUTPUT=""
TEMP=""
chain=""
USAGE="Usage: $0 -I <Input bam file> -g <global config> \"#:#-#\"]"
ERRORMESSAGE="#### ERROR"
ERRORMESSAGE1="The following error has occurred"

while getopts I:L:m:g:o:c:h:A: opt
  do     
  case "$opt" in
      I)      INP="$OPTARG";;
      L)      CHR="$OPTARG";;
      g)      GLOBAL="$OPTARG";;
      m)      MEM="$OPTARG";;
      o)      OUTPUT="$OPTARG";;
      c)      chain="$OPTARG";;
      A)      AUTO="$OPTARG";;
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


TEMP=$OUTPUT"_tempDir"

if [ ! -d $TEMP ]; then
    mkdir -p $TEMP
fi

JAVA="java -Xmx${HEAP}m -Djava.io.tmpdir="${TEMP}
GATK="$JAVA -jar "${GATKJAR}

if [[ $REFTYPE == "hg" ]]  # hg18/19 -> chr1, chr2 etc;  build36/37 -> 1, 2 etc                                            
    then
    cat $ExonFile | awk '{ print $1":"$2"-"$3}' > $TEMP"/target.list"
else
    cat $ExonFile |awk '{print $1":"$2"-"$3}' | sed 's/chr//' > $TEMP"/target.list"
fi

Targets=$TEMP"/target.list"

## dedup first
# use caution! 
# $SAMTOOLS rmdup $INP $INP.temp
# mv $INP.temp $INP
# $SAMTOOLS index $INP


$GATK \
    -R $REF \
    --DBSNP $DBSNP \
    -I $INP \
    -L $Targets \
    -nt 2 \
    -T CountCovariates \
    -cov ReadGroupCovariate \
    -cov QualityScoreCovariate \
    -cov CycleCovariate \
    -cov DinucCovariate \
    -recalFile $INP.recal_data.csv

# [-B:mask,VCF sitesToMask.vcf] \

if [[ $OUTPUT == ""  ]]
    then
    OUTPUT=$INP.recalibrated.bam
fi

$GATK \
    -R $REF \
    -T TableRecalibration \
    -I $INP \
    -recalFile $INP.recal_data.csv \
    -compress 5 \
    --out $OUTPUT

### need to 
## according to lh3 at http://biostar.stackexchange.com/questions/1269/what-is-the-best-pipeline-for-human-whole-exome-sequencing

# $SAMTOOLS calmd -br $OUTPUT.temp $REF > $OUTPUT
## update: this has been incorporated in GATK

$SAMTOOLS index $OUTPUT

rm -rf $TEMP
echo "recalibration complete"

d1=`dirname $INP`
job_name=`basename $d1`
if [[ $chain != ""  ]]; then
    # trigger downstream analysis 
    date
	if [[ $AUTO == "" ]]; then
	    cmd="qsub -N depth.$job_name -l mem=5G,time=30:: -o $OUTPUT.depth.log.o -e $OUTPUT.depth.log.e ${BPATH}/gatk_depthofcoverage.sh  -g $GLOBAL -I $OUTPUT  -m 3000 "
	else
	    cmd="qsub -N depth.$job_name.AUTO -l mem=5G,time=30:: -o $OUTPUT.depth.log.o -e $OUTPUT.depth.log.e ${BPATH}/gatk_depthofcoverage.sh  -g $GLOBAL -I $OUTPUT  -m 3000 -A AUTO"
	fi
    echo $cmd
    $cmd
   
fi

#
