#!/bin/bash
#$ -cwd

RunDir=$1
OutDir=$2

setting=$3

OLB=/ifs/data/c2b2/ngs_lab/ngs/usr/OLB-1.9.3/
PIPEBASE=/ifs/home/c2b2/ngs_lab/ngs/code/NGS/Pipeline/
status=/ifs/data/c2b2/ngs_lab/ngs/status/basecall.complete


if [[ $setting != "" ]]; then
	. $setting 
fi

BCALL=$RunDir/Data/Intensities/BaseCalls/

# THIS_QSEQS=$PROJECTS/$1/Qseqs
if [[ ! -e $OutDir ]]; then
    mkdir $OutDir
fi

chmod 750 $OutDir
pushd $OutDir

# cd $THIS_QSEQS

# check if the run is indeed finished
if [[ ! -s $RunDir/Basecalling_Netcopy_complete.txt || ! -s  $RunDir/RTAComplete.txt ]]; then 
    echo "Run not finished"
    exit 1
fi


absIN=`readlink -f $RunDir`
absOUT=`readlink -f $OutDir`
runName=`basename $absIN`

echo -e "Convert basecalls from $abspath to qseq file"


qseqout=$absOUT/qseq
fastqout=$absOUT/fastq

if [[ ! -e $fastqout ]]; then
    mkdir $fastqout
fi

if [[ ! -e $qseqout ]]; then
    mkdir $qseqout
fi

cd $qseqout

$OLB/bin/setupBclToQseq.py -b $BCALL -o $qseqout --overwrite -P .clocs 

make -j 8


if [[ -e "finished.txt" ]]; then
    echo -e "Qseq conversion from $absIN done" > BclToQseq.complete.txt
else
    echo "BclToQseq failed"
    exit 1
fi

echo -e "Convert qseq to fastq"

lanes=`ls s_*finished.txt`
for r in 1 2 3
  do 
  for f in $lanes
### example: ls s_*finished.txt
# s_1_finished.txt  s_3_finished.txt  s_5_finished.txt  s_7_finished.txt

    do
    lane=`echo $f | cut -f2 -d '_'`
    
    # run in background (a new process)
    cat s_"$lane"_"$r"_*_qseq.txt | $PIPEBASE/qseq2fastq.pl > $fastqout/s_"$lane"_"$r".fastq &   
  done
  wait  ##  wait for all background job to finish to start another loop 
  
done

touch $fastqout/QseqToFastq.complete.txt

echo -e "qseq to fastq done"

## clean up

bzip2 s*qseq.txt 

popd 
echo -e "conversion done" > $absOUT/BclToFastq.complete.txt
