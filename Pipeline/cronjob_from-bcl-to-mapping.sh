#!/bin/bash
#$ -cwd

setting="/ifs/home/c2b2/ngs_lab/ngs/code/NGS/Pipeline/global_setting.sh"

. $setting

#HiSeqRuns=/ifs/scratch/c2b2/ngs_lab/ngs/Illumina/
#FastqDir=/ifs/scratch/c2b2/ngs_lab/ngs/Fastq/
#status=/ifs/data/c2b2/ngs_lab/ngs/status/basecall.complete
#PIPEBASE=/ifs/home/c2b2/ngs_lab/ngs/code/NGS/Pipeline/

runs=`ls $HiSeqRuns`

for f in $runs
  do
# echo $f
  
  g=`grep -w $f $StatusDir/basecall.complete`
  
  if [[ $g == "" ]]; then
      echo "process $f"
      mkdir -p $FastqDir/$f/
      # Bcl to fastq
      cmd="$PIPEBASE/bclToFastq.sh $HiSeqRuns/$f/ $FastqDir/$f/ $setting "
      $cmd
      echo $f >> $StatusDir/basecall.complete  
      echo "$cmd" >> $StatusDir/history.txt
      
      # demultiplex
      cmd="$PIPEBASE/demultiplex.sh $FastqDir/$f/fastq/ $FastqDir/$f/demultiplex/ $sampelSheet"
      $cmd
      echo "$cmd" >> $StatusDir/history.txt
  else
      echo "$f is processed"
  fi
  
done

