#!/bin/bash
#$ -cwd

if [[ $1 == "" ]];then
	echo "Missing argument Directory"
	exit
else
	DIR=$1
fi

EXOMEBASE="/ifs/data/c2b2/ngs_lab/ngs/code/NGS/Exome/"
UTILS="/ifs/data/c2b2/ngs_lab/ngs/code/NGS/utils/"


# summary stats on depth: 
## reads mapped to whole genome
cd $DIR
if [ -e $DIR/release ];then
	suff=`date +%N`
	release="$DIR/release_$suff"
else
	release="$DIR/release"
fi
mkdir $release
cd $release

#SNS release (PASS filter)
sh $EXOMEBASE/do_release.sh $DIR/VarCalling/list.vcf-files.txt.vcf.complete.annotated.vcf.filtered.vcf.snv
# indel release:
sh $EXOMEBASE/do_release.sh $DIR/VarCalling/list.vcf-files.txt.vcf.complete.annotated.vcf.indel.filtered.vcf
ln -s $DIR/VarCalling/list.vcf-files.txt.vcf.complete.annotated.vcf raw.variants.vcf
ln -s $DIR/VarCalling/list.vcf-files.txt.vcf.complete.annotated.vcf.filtered.vcf.snv.release SNVs.vcf
ln -s $DIR/VarCalling/list.vcf-files.txt.vcf.complete.annotated.vcf.indel.filtered.vcf.release indels.vcf 

# Polyphen2 and SIFT scores. (use annovar)
qsub -o $DIR/logs/annovar.o -e $DIR/logs/annovar.e $EXOMEBASE/do_annovar.sh $DIR/$release/SNV.vcf


## SNV summary stats:
# coding regions
ruby $UTILS/vcf_transition-transversion-per-sample.rb  SNVs.vcf > summary.coding.xls 
# non-coding regions
ruby $UTILS/vcf_transition-transversion-per-sample.rb  SNVs.vcf -1 > summary.noncoding.xls 

for f in `ls $DIR/mapping/*refine/all.recalibrated.bam`; do echo $f; done > list.recalib_bam.list.release

echo "" >summary.a
echo "" >summary.b
echo "" >summary.c
echo "" >summary.d

for f in `cat list.recalib_bam.list.release`; do 
	head -2 $f.coverage.sample_summary | tail -1 |cut -f1 >>summary.a
	mappedreads=`head -1 $f.reads.mapped | sed 's/ \+/\t/'|cut -f1`
	echo $mappedreads >>summary.b
	head -2 $f.coverage.sample_summary | tail -1 |cut -f2,3,5,7,8,9,10,11  >>summary.c
	targetreads=`head -2 $f.coverage.sample_summary | tail -1 |cut -f2 ` 
	perc=$(($targetreads/$mappedreads))
	echo $perc >> summary.d
done

## information about the targeted regions.
echo -e "sample_id\t#Raw_Reads\t#Bases_mapped_on_Target\tAvg_depth\tMedian_depth\tD1(%)\tD5(%)\tD10(%)\tD15(%)\tD20(%)\t%Reads_mapped_on_target" > summary_reads.xls
paste summary.a summary.b summary.c summary.d >> summary_reads.xls

rm summary.a
rm summary.b
rm summary.c
rm summary.d


tem=`dirname $DIR`
project=`basename $tem`
samp=`wc -l list.recalib_bam.list.release | awk '{print $1;}'`

echo -e "This data is released at http://genomecenter.columbia.edu/ngs/  \n" > README
echo -e "Date : `date +%x` \n" >> README
echo -e "ProjectID : $project \n" >> README
echo -e "Note: \n" >> README
echo -e "   - No. of Samples = 12 " >> README
echo -e "   - There are 6 files in all. " >> README
echo -e "   - raw.variants.vcf  : Contains all raw variants " >> README
echo -e "   - indels.vcf 	: Contains filtered indels  " >> README
echo -e "   - SNVs.vcf    	: Contains filtered SNV calls" >> README
echo -e "   - summary.coding	: Summary files " >> README
echo -e "   - summary.noncoding	: Summary files " >> README
echo -e "   - summary_reads	: Summary of reads, coverage, D15 etc. per sample \n" >> README

echo -e "  Data on the web will be available for a month, after that, the data must be requested again." >> README

echo -e "  Pease write to sequencing@columbia.edu for questions. \n" >> README

rm  list.recalib_bam.list.release

