#!/bin/bash
#$ -cwd
# Findmem

# Batch Query by chromosme for EVS 

EVSJAR="/ifs/scratch/c2b2/ngs_lab/sz2317/softwares/evsClient-v.0.0.2/evsClient.jar"

HEAP=1
CHR=""
BEGIN=""
END=""
OUTDIR=""
USAGE="Usage: $0 outdirectory chromosome begin end "
ERRORMESSAGE="#### ERROR"

OUTDIR=$1
# [CHR, BEGIN, END] = split($2)
arr=($(echo $2 | tr " " "\n"))

CHR=${arr[0]}
BEGIN=${arr[1]}
END=${arr[2]}


if [[ $CHR == "" ]] || [[ $END == "" ]] || [[ $BEGIN == "" ]]; then
	echo $USAGE
	exit 1
fi

if [[ $OUTDIR == "" ]]; then
	$OUTDIR=`pwd`
#echo $USAGE
	OUTDIR=$OUTDIR"/"
#	exit 1 
fi

if [ ! -d $OUTDIR ]; then
    mkdir $OUTDIR
fi

TEMP=$OUTDIR"/temp/"

if [ ! -d $TEMP ]; then
    mkdir -p $TEMP
fi

JAVA="java -Xmx${HEAP}g -Djava.io.tmpdir="${TEMP}
EVSBATCH="$JAVA -jar "${EVSJAR}

#if [[ $CHR == "" ]]
#    then
#    echo "no chr prodived"
#fi
#if [[ $CHR == "23" ]]
#    then
#    CHR="X"
#fi
#if [[ $CHR == "24"  ]]
#    then
#    CHR="Y"
#fi

#echo $CHR 

#1       14467	228711
$EVSBATCH  -t $CHR:$BEGIN-$END -f vcf

rm -fr $TEMP
