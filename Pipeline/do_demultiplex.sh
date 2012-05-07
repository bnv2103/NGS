#!/bin/sh
#$ -l mem=8G,time=8::

RUBY18="/ifs/data/c2b2/ngs_lab/ngs/usr/local/bin/ruby"
NGSSHELL="/ifs/data/c2b2/ngs_lab/ngs/code/shell/"
PIPEBASE="/ifs/data/c2b2/ngs_lab/ngs/code/NGS/Pipeline/"

fastqout=$1
demultiplexout=$2
outprefix=$3
runName=$4

cmd="$RUBY18 $PIPEBASE/demultiplex.rb $fastqout $demultiplexout $outprefix $demultiplexout/$runName.tsv 8"
$cmd

echo "Cluster Density PF for all lanes"
absOUT=`dirname $fastqout`
awk 'BEGIN{lane=0;ct=0;sum=0} /^[0-9]/{if ($1==lane) {ct++; sum+=$3;} else{ print lane, sum/ct; lane=$1; ct=1;sum=$3;}}END{ print lane, sum/ct;}'  $absOUT/reports/NumClusters\ By\ Lane\ PF.txt 

## clean up
##bzip2 s*qseq.txt  & 
echo -e "conversion done" > $absOUT/BclToFastq.complete.txt

cd $demultiplexout

#Trigger Automatic Pipeline
cmd="sh $PIPEBASE/post_demux.sh $demultiplexout $runName "
echo $cmd
$cmd
	

