#!/bin/bash
#$ -cwd
#$ -l mem=6G,time=8::

#Run as  
#  qsub -o /ifs/scratch/c2b2/ngs_lab/ngs/Projects/WG-seq/120309_ANDREA_GABRIELLE_6_HUMAN_GENOME_30X_PE_HISEQ/120316_SN650_0258_AC0EAMACXX/logs/split.8_3.o -e /ifs/scratch/c2b2/ngs_lab/ngs/Projects/WG-seq/120309_ANDREA_GABRIELLE_6_HUMAN_GENOME_30X_PE_HISEQ/120316_SN650_0258_AC0EAMACXX/logs/split.8_3.e ../split_fastq.sh   s_8_3.fastq s_8_2.total AC1 120316_SN650_0258_AC0EAMACXX

# Input ex. s_2_1.fastq.bz2 OR unzipped 
# output prefix example  120224_SN650_0256_BD0PA5ACXX_lane7_71_1.fastq

lane=`basename $1 | cut -f2 -d"_" `
suf=`basename $1 | cut -f3 -d"_" | sed 's/.bz2//'`
pref=$4"_lane"	#RunName "120210_SN650_0253_AD0JGWACXX"
sample=$3
fq=`readlink -f $1`

barcode_file=`echo $fq | sed 's/_1.fast.\+/_2.fastq.barcode-stats/g' `
# lines=`tail -1 $2`	#TODO possibly change 2nd argument to no. of splits needed
tot=`cut -f2  $barcode_file | awk 'BEGIN{s=0}{s=s+$1;}END{print s;}' `
lines=`echo $(($tot/2/4*4)) `	#default split into 8 pieces.

echo "file $fq"
echo "lane -- $lane"
echo "suff  $suf"
echo "prefix  $pref"
echo "samples  $sample"
echo "lines $lines"
dir=$pref$lane"_"$sample"_"$suf
echo "dir $dir"

if [[ ! -e $dir ]]; then
	mkdir $dir
fi
cd $dir

if [[ $fq =~ ".bz2" ]]; then
	bzcat $fq | split --lines=$lines &
else
	split --lines=$lines $fq &
fi

cd ..

#If Paired End then split those files.
fq3=`echo $fq | sed 's/_1.fastq/_3.fastq/' `
if [ -e $fq3 ]
then
	suf=`basename $fq3 | cut -f3 -d"_" | sed 's/.bz2//'`
	dir3=$pref$lane"_"$sample"_"$suf
	echo "dir $dir3"
	echo "file $fq3"
	if [[ ! -e $dir3 ]]; then
	        mkdir $dir3
	fi
	cd $dir3
	if [[ $fq3 =~ ".bz2" ]]; then
	        bzcat $fq3 | split --lines=$lines
	else
        	split --lines=$lines $fq3
	fi
fi

cd ..
echo "split complete"
wait

#Fire Mapping
sh /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/do_mapping.sh $dir /ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/global_setting_b37.sh
