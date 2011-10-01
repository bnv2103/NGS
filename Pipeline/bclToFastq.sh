#!/bin/bash
#$ -cwd


OLB=/ifs/data/c2b2/ngs_lab/ngs/usr/OLB-1.9.3/
PIPEBASE=/ifs/home/c2b2/ngs_lab/ngs/code/NGS/Pipeline/
StatusDir=/ifs/data/c2b2/ngs_lab/ngs/status/
RUBY18="/ifs/data/c2b2/ngs_lab/ngs/usr/local/bin/ruby"
setting="/ifs/home/c2b2/ngs_lab/ngs/code/NGS/Pipeline/global_setting.sh"


USAGE="Usage: $0 -i RunDir  -o OutDir  -s setting [-n num_threads]"

nt=8 # default 8 threads for bcl to qseq

while getopts i:o:s:n:h opt
  do
  case "$opt" in
      i) RunDir="$OPTARG";;
      o) OutDir="$OPTARG";;
      s) setting="$OPTARG";;
      n) nt="$OPTARG";;
      h) echo $USAGE; exit 1
  esac
done

if [[ $RunDir == "" || $OutDir == ""  || $setting == "" ]]
    then
    echo $USAGE
    exit 1
fi

. $setting 

BCALL=$RunDir/Data/Intensities/BaseCalls/
absIN=`readlink -f $RunDir`  ## not robust enough, need a better method
absOUT=`readlink -f $OutDir`
runName=`basename $absIN`

echo "outdir: $OutDir     $absOUT "
# echo $absIN

if [[ ! -e $absOUT ]]; then
    mkdir $absOUT
fi

chmod 750 $absOUT
pushd $absOUT

# check if the run is indeed finished
if [[ ! -s $absIN/Basecalling_Netcopy_complete.txt || ! -s  $absIN/RTAComplete.txt ]]; then 
    echo "Run not finished"
    exit 1
fi

echo -e "Convert basecalls from $absIN to qseq file"


qseqout=$absOUT/qseq
fastqout=$absOUT/fastq
demultiplexout=$absOUT/demultiplex

if [[ ! -e $fastqout ]]; then
    mkdir $fastqout
fi

if [[ ! -e $qseqout ]]; then
    mkdir $qseqout
fi

cd $qseqout

$OLB/bin/setupBclToQseq.py -b $BCALL -o $qseqout --overwrite -P .clocs 

/ifs/data/c2b2/ngs_lab/ngs/usr/bin/make  -j $nt
## alternatively: qsub  -pe make 8-10 qmake.sh 
# cat qmake.sh
# #!/bin/bash
# #$ -cwd

# qmake -inherit --

if [[ -e "finished.txt" ]]; then
    echo -e "Qseq conversion from $absIN done" > BclToQseq.complete.txt
else
    echo "BclToQseq failed"
    exit 1
fi

echo -e "Convert qseq to fastq"

### example: ls s_*finished.txt
# s_1_finished.txt  s_3_finished.txt  s_5_finished.txt  s_7_finished.txt
lanes=`ls s_*finished.txt`

reads="1 2"  # default SE
if [[ -e  $absIN/Basecalling_Netcopy_complete_Read3.txt  ]]; then
    reads="1 2 3"  # PE
fi

njobs=0
for r in $reads
  do 
  for f in $lanes
    do
    lane=`echo $f | cut -f2 -d '_'`
    
    # run in background (a new process)
    $PIPEBASE/qseq2fastq.pl s_"$lane"_"$r"_*_qseq.txt > $fastqout/s_"$lane"_"$r".fastq &   
    let njobs=$njobs+1

    if [[ $njobs == $nt ]]; then
	wait  ##  wait for all background job to finish to start another loop 
	njobs=0
    fi
  done
done

wait


touch $fastqout/QseqToFastq.complete.txt

echo -e "qseq to fastq done"

## do  demultiplexing

# demultiplex                                                                                                     
sampleSheet=$SampleSheets/$runName.csv

if [[ -s $sampleSheet ]]; then
    

    if [[ ! -e $demultiplexout ]]; then
	mkdir $demultiplexout
    fi
    cp $sampleSheet $demultiplexout/$runName.csv
    outprefix=`echo $runName | cut -f1 -d '_' `
    cmd="$RUBY18 $PIPEBASE/demultiplex.rb $fastqout $demultiplexout $outprefix $demultiplexout/$runName.csv $nt"
    $cmd
    echo "$cmd" >> $StatusDir/history.txt
else
    echo "SampleSheet is missing. Failed to start demultiplexing."
fi

## clean up
## bzip2 s*qseq.txt  & 


popd 
echo -e "conversion done" > $absOUT/BclToFastq.complete.txt



