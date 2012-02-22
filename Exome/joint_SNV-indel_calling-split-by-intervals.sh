#!/bin/bash
#$ -cwd

nt=1 # default number of threads
dcov=300 # down sampling to dcov if depth is larger than dcov
njobs=100
USAGE="Usage: $0 -i <list of bam files> -m <heap> -s <global setting> [ -n number_of_threads] [ -j number_of_qjobs] [-d down_sampling] [ -v total_mem ] [ -o output_Dir ]"

while getopts i:m:o:s:d:n:j:b:q:v:h:t:A: opt
  do 
  case "$opt" in
      i) bamlist="$OPTARG";;
      m) MEM="$OPTARG";;
      o) outDir="$OPTARG";;
      s) setting="$OPTARG";;  # global config
      n) nthreads="$OPTARG";;
      d) dcov="$OPTARG";;
      j) njobs="$OPTARG";;
      t) qtime="$OPTARG";;
#      b) mbq="$OPTARG";;
#      q) mmq="$OPTARG";;
      v) qmem="$OPTARG";;
      A) AUTO="$OPTARG";;
      h) echo $USAGE
	  exit 1;;
  esac
done

if [[ $bamlist == ""  || $setting == "" ]]
    then
    echo $USAGE
    exit 1
fi

d1=`dirname $bamlist`
basen=`basename d1`
if [[ $AUTO == ""  ]];then
	job_ext=".$basen"
else
	job_ext=".$basen.AUTO"
fi

if [[ ! $nthreads == "" ]]; then
    nt=$nthreads
fi

let heap=$nt*4

if [[ ! $MEM == "" ]]; then
    heap=$MEM
fi

if [[ $qmem == "" ]]; then
    let qmem=$heap+3  # allocated RAM for each qjob in total
fi

if [[ $qtime == "" ]]; then
    let qtime=240  # time limit for the job
fi

. $setting

### get the name convention for chromosomes
# chrString=`head -2 $REF | egrep "^>chr"`
#if [[ $chrString != "" ]]; then
#    chrprefix='chr'
#else
#    chrprefix=''
# fi

bamlist=`readlink -e $bamlist`

temp=$bamlist"_VarCalling_dir"

if [[ $outDir != "" ]]; then
	temp=$outDir
fi

mkdir -p $temp

# echo $ExonFile

target=$temp"/all-targets.list"

awk '{print $1":"$2"-"$3}' $ExonFile > $target
    
num=`wc -l $target | awk '{print \$1}'`

let nslice=$num/$njobs+1

bname=`basename $bamlist`
total=0
for (( j=1; j<=$njobs; j++ ))  #  
  do
  tempd=${temp}"/slices/"$j

  mkdir -p $tempd
  
  out=${temp}"/var.slice."$j".sh"
#  chromosome=${chrprefix}${i}
  chrtarget=$out".targets.list" 

  let total=$total+$nslice

  if [[ $j -eq $njobs ]]
      then
      let njobs=$njobs-1
      let nslice=$num-$nslice*$njobs
      let njobs=$njobs+1
  fi
 
  head -${total} $target | tail -${nslice} > $chrtarget

  infofields="-A AlleleBalance  -A DepthOfCoverage  -A BaseQualityRankSumTest  -A HomopolymerRun -A MappingQualityRankSumTest -A MappingQualityZero -A QualByDepth  -A RMSMappingQuality  -A SpanningDeletions  -A HaplotypeScore "
  
  echo '#!/bin/bash'  > $out
  echo '#$ -cwd' >> $out
  echo 'uname -a' >> $out
  cmd="java -Xmx${heap}g -Djava.io.tmpdir=${tempd}  -jar $GATKJAR -T UnifiedGenotyper  -R $REF   -nt ${nt} -o ${temp}/var.slice.$j.raw.vcf -stand_call_conf 50.0 -stand_emit_conf 10.0 -dcov ${dcov} -glm BOTH  -L $chrtarget -I $bamlist -metrics ${temp}/var.slice.$j.raw.vcf.metrics -G Standard  -B:dbsnp,VCF ${DBSNPVCF} -B:compdbSNP132,VCF $DBSNP132 $infofields"
  
  echo $cmd >> $out
  echo "if [ ! -e $temp/status.Varcalling ];then touch $temp/status.Varcalling; fi " >> $out
  echo "echo $j >> $temp/status.Varcalling " >> $out
  echo "completed=\`wc -l $temp/status.Varcalling | awk '{print \$1}'\` " >> $out
  echo "while [[ \$completed != \"\" ]]; do " >> $out
  echo "if [[ \$completed -eq \"100\" ]];then  echo \"all completed\" >>  $temp/status.Varcalling; break;" >> $out
  echo "elif [[ \$completed -lt \"95\" || \$completed -gt \"100\" ]]; then exit; fi" >>$out
  echo "else " >> $out
  echo "sleep 60 " >> $out
  echo "completed=\`wc -l $temp/status.Varcalling | awk '{print \$1}'\` " >> $out
  echo "fi; done" >> $out
  echo " for (( i=1;i <= 100;i++ )); do ls $outDir/var.slice.\$i.raw.vcf; done > $outDir/list.vcf-files.txt " >> $out
  echo " sh ${BPATH}/vcf_concat_slices.sh $outDir/list.vcf-files.txt $outDir/list.vcf-files.txt.vcf " >> $out

if [[ $AUTO == "" ]]; then  
  echo " qsub -l mem=6G,time=$qtime:: -N annotation$job_ext -o annotation.o -e annotation.e ${BPATH}/gatk_annotator.sh  -v $outDir/list.vcf-files.txt.vcf -g $setting -m 4 -b $bamlist  " >> $out
else
  echo " qsub -l mem=6G,time=$qtime:: -N annotation$job_ext -o annotation.o -e annotation.e ${BPATH}/gatk_annotator.sh  -v $outDir/list.vcf-files.txt.vcf -g $setting -m 4 -b $bamlist -A AUTO " >> $out
fi

  qsub -l mem=${qmem}G,time=${qtime}:: -o $temp/log.$j.o -e $temp/log.$j.e -N var.$j$job_ext $out 
  # echo $qmem
done
