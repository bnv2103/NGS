#!/bin/bash
#$ -cwd
#$ -l mem=6G,time=12::

#Run as  
#  qsub -o /ifs/scratch/c2b2/ngs_lab/ngs/Projects/WG-seq/120309_ANDREA_GABRIELLE_6_HUMAN_GENOME_30X_PE_HISEQ/120316_SN650_0258_AC0EAMACXX/logs/split.8_3.o -e /ifs/scratch/c2b2/ngs_lab/ngs/Projects/WG-seq/120309_ANDREA_GABRIELLE_6_HUMAN_GENOME_30X_PE_HISEQ/120316_SN650_0258_AC0EAMACXX/logs/split.8_3.e ../split_fastq.sh   s_8_3.fastq s_8_2.total AC1 120316_SN650_0258_AC0EAMACXX

# Input ex. s_2_1.fastq.bz2 OR unzipped 
# output prefix example  120224_SN650_0256_BD0PA5ACXX_lane7_71_1.fastq

BPATH="/ifs/scratch/c2b2/ngs_lab/sz2317/scripts/WGS/"
lane=`basename $1 | cut -f2 -d"_" `
suf=`basename $1 | cut -f3 -d"_" | sed 's/.bz2//'`
pref=$4"_lane"	#RunName "120210_SN650_0253_AD0JGWACXX"
sample=$3
fq=`readlink -f $1`

barcode_file=`echo $fq | sed 's/_1.fast.\+/_2.fastq.barcode-stats/g' `
splits=$2	#TODO possibly change 2nd argument to no. of splits needed , if 0 then just unzip without splitting, but keep same directory structure, 2=means split into 8 pieces, 1=means split into 4 pieces.
tot=`cut -f2  $barcode_file | awk 'BEGIN{s=0}{s=s+$1;}END{print s;}' `

workingdir=`readlink -f  $1 | dirname `
echo "Entering Working Directory $workingdir "

cd $workingdir

echo "file $fq"
echo "lane -- $lane"
echo "suff  $suf"
echo "prefix  $pref"
echo "samples  $sample"
dir=$pref$lane"_"$sample"_"$suf
echo "dir $dir"

if [[ ! -e $dir ]]; then
	mkdir $dir
fi
cd $dir

if [[ $splits == 0 ]]
then	## DOnot split 
	echo "Zero splits"
	if [[ $fq =~ ".bz2" ]]; then
        	bunzip2 -ck $fq > xaa &
	else
        	cp $fq . &
	fi
else
	lines=`echo $(($tot/$splits/4*4)) `	# if 0 then just unzip without splitting, but keep same directory structure, 2=means split into 8 pieces, 1=means split into 4 pieces.
	echo "lines $lines"

	if [[ $fq =~ ".bz2" ]]; then
		bzcat $fq | split --lines=$lines &
	else
		split --lines=$lines $fq &
	fi
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
	if [[ $splits == 0 ]]
	then    ## DOnot split
        	echo "Zero splits"
	        if [[ $fq3 =~ ".bz2" ]]; then
	                bunzip2 -ck $fq3 > xaa &
	        else	
	                cp $fq3 . &
        	fi
	else
		if [[ $fq3 =~ ".bz2" ]]; then
		        bzcat $fq3 | split --lines=$lines
		else
        		split --lines=$lines $fq3
		fi
	fi
fi

cd ..
echo "split complete"
wait

#Fire Mapping
sh $BPATH/do_mapping.sh $dir $BPATH/global_setting_b37.sh
