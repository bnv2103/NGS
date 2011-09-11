#!/bin/bash
#$ -cwd

###### To consider ######
# right now the status is recorded in a plain text file. We should
# consider to replace it with a sqlite database. 


setting="/ifs/home/c2b2/ngs_lab/ngs/code/NGS/Pipeline/global_setting.sh"

#HiSeqRuns=/ifs/scratch/c2b2/ngs_lab/ngs/Illumina/
#FastqDir=/ifs/scratch/c2b2/ngs_lab/ngs/Fastq/
#StatusDir=/ifs/data/c2b2/ngs_lab/ngs/status/basecall.complete
#PIPEBASE=/ifs/home/c2b2/ngs_lab/ngs/code/NGS/Pipeline/
#SampleSheets=/ifs/data/c2b2/ngs_lab/ngs/status/SampleSheets/
#QueueLog=/ifs/data/c2b2/ngs_lab/ngs/status/QueueLog/
. $setting


runs=`ls $HiSeqRuns`

for f in $runs
  do
# echo $f
  
  g=`grep -w $f $StatusDir/basecall.complete`
  
  if [[ $g == "" ]]; then
      # check if the run is complete
      if [[ -s $HiSeqRuns/$f/Basecalling_Netcopy_complete.txt && -s  $HiSeqRuns/$f/RTAComplete.txt ]]; then
	  echo "process $f"
	  echo -e "$f\t$HiSeqRuns/$f/\t$FastqDir/$f" >> $StatusDir/basecall.complete  

      # Bcl to fastq, will take 4-5 hours. then do demultiplexing
	  cmd="qsub -N process.$f -l mem=8G,time=72::  -pe smp 4   -o $QueueLog/bclToFastq.$f.o -e $QueueLog/bclToFastq.$f.e $PIPEBASE/bclToFastq.sh -i $HiSeqRuns/$f -o $FastqDir/$f -s $setting -n 4"
	  $cmd
	  echo $cmd
	  echo "$cmd" >> $StatusDir/history.txt
#	  popd
      else
	  echo "Run not finished"
      fi
  else
      echo "$f is processed"
  fi
done

