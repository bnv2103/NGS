#!/bin/bash
#$ -cwd
# Findmem


HEAP=7000


INP=""
CHR=""
List=""
ExonFile=""
TEMP=""
MEM=""

USAGE="Usage: $0 -I <Input bam file> -R <Reference fasta> [-L \"#:#-#\"] [-E <ExonFile>]"
EXONUSAGE="Please specify either the file containing the interval list using -E or the sequences using -L"

while getopts I:L:g:t:m:h:A: o
do      case "$o" in
        I)      INP="$OPTARG";;
        g)      GLOBAL="$OPTARG";;
        t)      TEMP="$OPTARG";;
        m)      MEM="$OPTARG";;
	A)	AUTO="$OPTARG";;
        h)      echo $USAGE
                exit 1;;
        esac
done

if [[ $INP == "" || $GLOBAL == ""  ]]
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
if [ $JOB_ID == "" ]; then
    JOB_ID="depth"
fi

if [[ $TEMP == "" ]]; then
    TEMP=$INP"_depth_temp"
fi

if [ ! -d $TEMP ]; 
    then    
    mkdir -p $TEMP
fi

JAVA="java -Xmx${HEAP}m -Djava.io.tmpdir="${TEMP}
GATK="$JAVA -jar "${GATKJAR}

CHR="$INP.target.list"

    
if [[ $REFTYPE == "hg" ]]  # hg18/19 -> chr1, chr2 etc;  build36/37 -> 1, 2 etc                                            
    then
    cat $ExonFile | awk '{ print $1":"$2"-"$3}' > "$INP.target.list"
else
    cat $ExonFile | awk '{print $1":"$2"-"$3}' | sed 's/chr//' > "$INP.target.list"
fi
	
$GATK \
 -T DepthOfCoverage \
 -L $CHR \
 -R $REF \
 -I $INP \
 -o $INP.coverage \
 -ct 1 \
 -ct 5 \
 -ct 10 \
 -ct 15 \
 -ct 20

rm -rf $TEMP

if [ $? == 0 ]
    then
    rm -f $INP.coverage
    echo "depth of coverage complete"
fi

# Could possibly include these

# Summary stats on depth:
## reads mapped to whole genome
echo $INP > $INP.reads.mapped
samtools flagstat $INP  >> $INP.reads.mapped
samtools idxstats $INP > $INP.idxstats
## information about the targeted regions.
echo -e "sample_id\ttotal\tmean\tgranular_third_quartile\tgranular_median\tgranular_first_quartile\tD1\tD5\tD10\tD15\tD20" > $INP.reads.target
cat $INP.coverage.sample_summary  | grep -v total | grep -v Total >> $INP.reads.target 

if [[ $AUTO != "" ]]
then
        #Trigger automatic downstream steps : joint var calling
        readlink -f $INP > $INP.list
	path_inp=`dirname $INP`
       cmd_varcalling="sh ${BPATH}/joint_SNV-indel_calling-split-by-intervals.sh -i $INP.list -m 8 -s $GLOBAL -n 2 -j 100 -d 300 -v 10 -t 30 -o $path_inp/VarCalling -A AUTO "
	echo $cmd_varcalling
	$cmd_varcalling	
fi

