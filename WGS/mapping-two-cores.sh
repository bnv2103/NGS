#!/bin/bash
#$ -cwd

## NOTE: this require two cores. For titan as of June 2011, to get two cores for bwa requires set "-l mem=8G" or larger. 

# default values
maxgaps=2
maxeditdist=0.04
qualtrim=5
platform="illumina"
threads=2
sampleName=""
readgroup=""
bwa=`which bwa`
samtools=`which samtools`
output=""
bwaversion=`$bwa 2>&1 | grep Version | awk '{print $2""$3}'`
#sortmem=1000000000  # mem allocated for samtools sort  
chain="0"  # whether call downstream analysis (default 0 means not). 
setting=""

### Note: bwa sample uses about 3.5G RAM

USAGE="Usage: $0 -i foo_1.fastq  -s global_setting [ -p foo_2.fastq ] [ -g maxgaps] [ -q qualtrim ] [ -z readgroup] [ -n sampleName] [ -f platform] [-o output_prefix]"

while getopts i:p:g:q:d:n:t:s:z:f:m:o:c:h:A:y: opt
  do      
  case "$opt" in
      i) fastq1="$OPTARG";;
      p) fastq2="$OPTARG";;
      m) sortmem="$OPTARG";;
      g) maxgaps="$OPTARG";;
      d) maxeditdist="$OPTARG";;
      q) qualtrim="$OPTARG";;
      n) sampleName="$OPTARG";;
      t) threads="$OPTARG";;
      z) readgroup="$OPTARG";;
      y) ID="$OPTARG";;
      f) platform="$OPTARG";;
      o) output="$OPTARG";;
      c) chain="$OPTARG";;
      s) setting="$OPTARG";;
      A) AUTO="$OPTARG";;
      h)    echo $USAGE
          exit 1;;
  esac
done

if [[ $fastq1 == "" || $setting == "" ]]; then
    echo $USAGE
    exit 1
fi

date

. $setting

if [[ $readgroup == "" ]]; then
    readgroup=`basename $fastq1 | sed 's/.fastq$//'  | sed 's/.txt$//'`
fi

if [[ $sampleName == "" ]]; then
    sampleName=$readgroup
fi

if [[ $ID == "" ]]; then
    ID=$readgroup
fi

if [[ $output == "" ]]; then
    output=$fastq1.sorted
fi


## read group specification:
##          -r STR   read group header line such as `@RG\tID:foo\tSM:bar' [null]
rgheader="@RG\tID:$ID\tSM:$sampleName\tLB:$readgroup\tPL:$platform\tCN:NGSColumbia"

######### align step
cmd="$bwa aln -I -q $qualtrim -o $maxgaps -n $maxeditdist -t  $threads  $REF  $fastq1 > $fastq1.sai"
echo $cmd
$bwa aln -I -q $qualtrim -o $maxgaps -n $maxeditdist -t  $threads  $REF  $fastq1 > $fastq1.sai  &  ### run in background


if [[ ! $fastq2 == "" ]]; then  # paired-ends
    cmd="$bwa aln -I -q $qualtrim -o $maxgaps  -n $maxeditdist -t  $threads  $REF  $fastq2 > $fastq2.sai"
    echo $cmd
    $bwa aln -I -q $qualtrim -o $maxgaps  -n $maxeditdist  -t  $threads  $REF  $fastq2 > $fastq2.sai 

# view -bS -o #{output} -

   # $bwa sampe -p $platform -i $readgroup -l $readgroup -m $sampleName $ref $fastq1.sai $fastq2.sai $fastq1 $fastq2 | $samtools view -bS - | $samtools sort -m $sortmem  -  $output

    
    wait   ### need to wait the F reads finish

    date
    cmd="$bwa sampe -P -r $rgheader  $REF $fastq1.sai $fastq2.sai $fastq1 $fastq2 | $samtools view -bS -  > $output.bam.temp"
    echo $cmd
    $bwa sampe -P -r $rgheader  $REF $fastq1.sai $fastq2.sai $fastq1 $fastq2 | $samtools view -bS -  > $output.bam.temp

    date
    echo "bwa alignment complete. Sorting the bam file ..."
    $samtools sort $output.bam.temp  $output
#    rm -f $fastq1.sai $fastq2.sai $output.bam.temp
   
else  # single-end
    
    wait
    
    date
    cmd="$bwa samse  -r $rgheader $REF $fastq1.sai $fastq1 | $samtools view -bS - >  $output.bam.temp"
    echo $cmd
    $bwa samse  -r $rgheader $REF $fastq1.sai $fastq1 | $samtools view -bS - >  $output.bam.temp
#    $samtools sort -m $sortmem  $output.bam.temp  $output
    date
    echo "bwa alignment complete. Sorting the bam file ..."
    $samtools sort $output.bam.temp  $output
   
#    rm -f $fastq1.sai $output.bam.temp
fi

# fix a bug in samtools sort / bwa
## $samtools view -H $output.bam | sed 's/SO\:unsorted/SO:coordinate/' > $output.bam.header

echo -e "@HD\tVN:1.0\tGO:none\tSO:coordinate" >  $output.bam.header
$samtools view -H $output.bam | egrep -v '^\@HD' >> $output.bam.header

# alternative for older version of echo:
# echo  "@PG"$'\t'"ID:$readgroup"$'\t'"VN:$bwaversion"$'\t'"CL:$bwa" >> $output.bam.header

$samtools reheader $output.bam.header $output.bam > $output.bam.temp

mv $output.bam.temp $output.bam
rm -f $output.bam.temp*

# index
$samtools index $output.bam
$samtools flagstat $output.bam > $output.flagstat
$samtools idxstats  $output.bam >  $output.idxstats

date
if [ -e $output.bam ] 
then
fqbase=`dirname $fastq1`
fqbase1=`dirname $fqbase`

zipcmd="qsub -o $fqbase1/zip.fq.o -e $fqbase1/zip.fq.e -l mem=1G,time=6:: $NGSSHELL/do_bzip2.sh  $fastq1  $fastq2 "
$zipcmd
echo $zipcmd
echo "splitting bam by chr"
sh $BPATH/splitChr_bam.sh  $output.bam &> $output.splitChr.log
fi
exit 


