#!/bin/bash
#$ -cwd

#e.g combined.complete.annotated.vcf

if [[ $1 == "" ]];then
	echo "Missing argument raw VCF file"
	exit
else
	vcf=`readlink -f $1`
	DIR=`dirname $vcf`
fi

EXOMEBASE="/ifs/data/c2b2/ngs_lab/ngs/code/NGS/Exome/"
UTILS="/ifs/data/c2b2/ngs_lab/ngs/code/NGS/utils/"
WebRelease="/ifs/data/shares/deepsequencing/solid/public_html/"
PIPEBASE="/ifs/data/c2b2/ngs_lab/ngs/code/NGS/Pipeline/"

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

 ln -s $vcf raw.variants.vcf
 ln -s $vcf.filtered.vcf.snv.release SNVs.vcf
 ln -s $vcf.indel.filtered.vcf.release indels.vcf 

## SNV summary stats:
# coding regions
 ruby $UTILS/vcf_transition-transversion-per-sample.rb  SNVs.vcf > summary.coding.xls 
# non-coding regions
 ruby $UTILS/vcf_transition-transversion-per-sample.rb  SNVs.vcf -1 > summary.noncoding.xls 

for f in `cat $DIR/vcf.list ` ; do
	g=`echo $f | sed 's/_refine\/.\+/_refine\/all.recalibrated.bam/'` ; 
	echo `readlink -f $g `;
done > list.recalib_bam.list.release

echo "" >summary.a
echo "" >summary.b
echo "" >summary.c
echo "" >summary.d

for f in `cat list.recalib_bam.list.release`; do 
	head -2 $f.coverage.sample_summary | tail -1 |cut -f1 >>summary.a
        g=`echo $f | sed 's/_refine\/all.recalibrated.bam/.flagstat/'` ;
	mappedreads=`grep total $g | sed 's/ \+/\t/'|cut -f1`
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

#E.g. /ifs/scratch/c2b2/ngs_lab/ngs/Projects/Exome-seq/120611_ALI_YIFU_20_HUMAN_EXOME_60X_PE_HISEQ/120629_SN828_0140_AD13G3ACXX/mapping_redo/Pulse-2.bam_refine/VarCalling/list.vcf-files.txt.vcf

project=` echo $DIR | tr "/" "\t" | cut -f9  `
samp=`wc -l list.recalib_bam.list.release | awk '{print $1;}'`

echo -e "This data is released at http://genomecenter.columbia.edu/ngs/$project  \n" > README
echo -e "Date : `date +%x` \n" >> README
echo -e "ProjectID : $project \n" >> README
echo -e "Note: \n" >> README
echo -e "   - No. of Samples = $samp " >> README
echo -e "   - There are 6 files in all. " >> README
echo -e "   - raw.variants.vcf  : Contains all raw variants " >> README
echo -e "   - indels.vcf 	: Contains filtered indels  " >> README
echo -e "   - SNVs.vcf    	: Contains filtered SNV calls" >> README
echo -e "   - summary.coding	: Summary of filtered coding variants " >> README
echo -e "   - summary.noncoding	: Summary of filtered non-coding variants " >> README
echo -e "   - summary_reads	: Summary of reads, coverage, D15 etc. per sample \n" >> README

echo -e "  Data on the web will be available for a month, after that, the data must be requested again." >> README

echo -e "  Pease write to sequencing@columbia.edu for questions. \n" >> README

rm  list.recalib_bam.list.release

if [ ! -d $WebRelease/$project ]
then
	mkdir $WebRelease/$project 
	cp * $WebRelease/$project/
	sh $PIPEBASE/make_release_index-html.sh $WebRelease/$project/
	echo "Web Release Created!"
else
        echo "Web Release NOT created!"
fi

