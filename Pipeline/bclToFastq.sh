#!/bin/bash
#$ -cwd


OLB=/ifs/data/c2b2/ngs_lab/ngs/usr/OLB-1.9.3/
PIPEBASE=/ifs/home/c2b2/ngs_lab/ngs/code/NGS/Pipeline/
StatusDir=/ifs/data/c2b2/ngs_lab/ngs/status/
RUBY18="/ifs/data/c2b2/ngs_lab/ngs/usr/local/bin/ruby"
setting="/ifs/home/c2b2/ngs_lab/ngs/code/NGS/Pipeline/global_setting.sh"
NGSSHELL="/ifs/data/c2b2/ngs_lab/ngs/code/shell/"

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

chmod 770 $absOUT
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

#if [[ -e "finished.txt" ]]; then
#    echo -e "Qseq conversion from $absIN done" > BclToQseq.complete.txt
#else
#    echo "BclToQseq failed"
#    exit 1
#fi

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
    $PIPEBASE/qseq2fastq.pl s_"$lane"_"$r"_*_qseq.txt > $fastqout/s_"$lane"_"$r".fastq  2> $fastqout/s_"$lane"n_"$r".fastq &   
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
#Clean up Qseq directory
flag=1
for i in `seq 1 8`
do
	for j in `seq 1 2`
	do
		if [ ! -e $fastqout/s_"$i"_"$j".fastq ] || [ ! -s  $fastqout/s_"$i"_"$j".fastq ]
		then flag=0
		fi
	done
done
if [ $flag -eq "1" ]
then
 	cp -r $absIN/Data/reports/ $absOUT/
	mkdir $absOUT/RunData
	cp -r $absIN/runParameters.xml $absOUT/RunData/
	cp -r $absIN/RunInfo.xml $absOUT/RunData/
	cp -r $absIN/Config $absOUT/RunData/
	cp -r $absIN/Data/Intensities/config.xml $absOUT/RunData/
	cp -r $absIN/Data/Intensities/RTAConfiguration.xml  $absOUT/RunData/
	cp -r $absIN/InterOp/  $absOUT/RunData/
	cp -r $absIN/Recipe/ $absOUT/RunData/
	rm s_*.txt
	rm -rf Matrix/
	rm -rf Phasing/
	rm -rf Plots/
	rm -rf SignalMeans/
fi

# do QC
# TODO
mkdir $fastqout/logs
for i in `seq 1 8`
do
	j=1
	echo ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/fastq_QCplot.rb $fastqout/s_$i"_$j.fastq" $fastqout/s_$i"_$j.pdf" | qsub -o $fastqout/logs/QC.$i.$j.o -e $fastqout/logs/QC.$i.$j.e -l mem=4G,time=8:: -N QC_lane.$i.$j
	j=3
	if [ -e $fastqout/s_$i"_$j.fastq" ];then
		echo ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/fastq_QCplot.rb $fastqout/s_$i"_$j.fastq" $fastqout/s_$i"_$j.pdf" | qsub -o $fastqout/logs/QC.$i.$j.o -e $fastqout/logs/QC.$i.$j.e -l mem=4G,time=8:: -N QC_lane.$i.$j
	fi
done

# get barcode stats
qsub -o $fastqout/logs/barcode-stats.o -e $fastqout/logs/barcode-stats.e -l mem=3G,time=48:: -N barcode_stats.$runName $PIPEBASE/do-barcode-stats.sh $fastqout/s_*_2.fastq

## do  demultiplexing
# demultiplex                                                                                                     
sampleSheet=$SampleSheets/$runName.csv

if [[ ! -s $sampleSheet ]]; then
    sampleSheet=$SampleSheets/$runName.tsv
fi 

if [[ -s $sampleSheet ]]; then
    
    
    if [[ ! -e $demultiplexout ]]; then
	mkdir $demultiplexout
	mkdir $demultiplexout"_n" 
    fi
    cp $sampleSheet $demultiplexout/$runName.csv
    cp $sampleSheet $demultiplexout"_n"/$runName.csv
    outprefix=`echo $runName`		# outprefix=`echo $runName | cut -f1 -d '_' `

#Demultiplex not passed filter reads to get stats only
    cmd="$RUBY18 $PIPEBASE/demultiplex_n.rb $fastqout "$demultiplexout"_n $outprefix "$demultiplexout"_n/$runName.csv 4 &"
    echo "$cmd" >> $StatusDir/history.txt
    $cmd

#Demultipex PF reads
    cmd="$RUBY18 $PIPEBASE/demultiplex.rb $fastqout $demultiplexout $outprefix $demultiplexout/$runName.csv $nt"
    echo "$cmd" >> $StatusDir/history.txt
    $cmd

echo "" > Reads.summary
for sm in `ls $demultiplexout"/*summary*" `; do
	head -2 $sm | awk '{ for (i=1; i<=NF; i++)  { a[NR,i] = $i; }}NF>p { p = NF; }END { for(j=1; j<=p; j++) { str=a[1,j]; for(i=2; i<=NR; i++){ str=str"\t"a[i,j];} print str; }}' >> Reads.summary
done
echo "" > Reads_n.summary
for sm in `ls $demultiplexout"_n/*summary*" `; do
	head -2 $sm | awk '{ for (i=1; i<=NF; i++)  { a[NR,i] = $i; }}NF>p { p = NF; }END { for(j=1; j<=p; j++) { str=a[1,j]; for(i=2; i<=NR; i++){ str=str"\t"a[i,j];} print str; }}' >> Reads_n.summary
done


else
    echo "SampleSheet is missing. Failed to start demultiplexing."
    exit
fi

echo "Cluster Density PF for all lanes"
awk 'BEGIN{lane=0;ct=0;sum=0} /^[0-9]/{if ($1==lane) {ct++; sum+=$3;} else{ print lane, sum/ct; lane=$1; ct=1;sum=$3;}}END{ print lane, sum/ct;}'  $absOUT/reports/NumClusters\ By\ Lane\ PF.txt 


## clean up
##bzip2 s*qseq.txt  & 

for i in `seq 1 8`; do 
 qsub -o $fastqout/zip."$i".o -e $fastqout/zip."$i".e -l mem=512M,time=4:: $NGSSHELL/do_bzip2.sh $fastqout/s_"$i"_1.fastq  $fastqout/s_"$i"n_1.fastq
 qsub -o $fastqout/zip."$i".o -e $fastqout/zip."$i".e -l mem=512M,time=2:: $NGSSHELL/do_bzip2.sh $fastqout/s_"$i"_2.fastq $fastqout/s_"$i"n_2.fastq
 if [[ -e $fastqout/s_"$i"_3.fastq ]];then
   qsub -o $fastqout/zip."$i".o -e $fastqout/zip."$i".e -l mem=512M,time=4:: $NGSSHELL/do_bzip2.sh $fastqout/s_"$i"_3.fastq  $fastqout/s_"$i"n_3.fastq
 fi
done  

#  popd 
echo -e "conversion done" > $absOUT/BclToFastq.complete.txt
cd $demultiplexout

#Trigger Automatic Pipeline
        cmd="sh $PIPEBASE/post_demux.sh $demultiplexout $runName "
        echo $cmd
        $cmd
	
	cd $demultiplexout

	touch "mailBody.txt"
        echo "" > "mailBody.txt"
        echo "Fastq Complete: The Fastq files for $f are ready ." >> mailBody.txt
        echo ""  >> mailBody.txt
        echo "Process: Automated DownStream Pipeline Begin " >> mailBody.txt
        echo "Hi-Seq Run: $f " >> mailBody.txt
        echo "Command: $cmd " >> mailBody.txt
        echo ""  >> mailBody.txt

 #       qstat -j process.$f >>  mailBody.txt
        cmd1="sh $PIPEBASE/sendMail.sh -t sz2317@c2b2.columbia.edu,xs2182@c2b2.columbia.edu,yshen@c2b2.columbia.edu,oc2121@c2b2.columbia.edu,ecb2152@c2b2.columbia.edu -s Hi-Seq-Run-$f-Complete -m mailBody.txt "
        echo $cmd1
        $cmd1
        rm mailBody.txt

popd
