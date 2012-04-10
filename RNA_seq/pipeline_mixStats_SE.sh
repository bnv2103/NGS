#!/bin/bash
#$ -cwd


#ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/printTitle.rb $fileName
#f1="111014_lane_2_48h-sample-2_1.fastq-se_cufflinks/accepted_hits.bam"
#f2="111014_lane_2_48h-sample-2_1.fastq-se_cufflinks_human/accepted_hits.bam"
f1=$1
f2=$2
#samtools view $f1 | cut -f1 | sort -u -S 16G > readsName1.txt
#samtools view $f2 | cut -f1 | sort -u -S 16G > readsName2.txt

#comm -12 readsName1.txt readsName2.txt > commReads.txt

samtools view -X $f1 | cut -f1,2 | grep 'pP' | cut -f1 |sort -u -S 16G > readsName1PE.txt
samtools view -X $f2 | cut -f1,2 | grep 'pP' | cut -f1 |sort -u -S 16G > readsName2PE.txt

# get overlaped reads
comm -12 readsName1PE.txt readsName2PE.txt > commReadsPE.txt
wc -l readsName1PE.txt > numOfSample1
wc -l readsName2PE.txt > numofSample2
wc -l commReadsPE.txt > numOfComm

#samtools view $f1 | cut -f1,5 > reads1.txt
#samtools view $f2 | cut -f1,5 > reads2.txt

#rm result1.txt
#rm result2.txt

#for i in `cat commReads.txt`; do grep -w $i reads1.txt >> result1.txt; done
#for i in `cat commReads.txt`; do grep -w $i reads2.txt >> result2.txt; done

#`Rscript compareMAPQ.R`










# ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/mix_comb_stats.rb $f1 $f2 $fileName
# ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/mpQC_test.rb $f1 $f2 $fileName

#f1="111014_lane_2_48h-sample-3_1.fastq-se_cufflinks"
#f2="111014_lane_2_48h-sample-3_1.fastq-se_cufflinks_human"

#ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/mix_comb_stats.rb $f1 $f2 $fileName

#f1="111014_lane_2_48h-sample-4_1.fastq-se_cufflinks"
#f2="111014_lane_2_48h-sample-4_1.fastq-se_cufflinks_human"

#ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/mix_comb_stats.rb $f1 $f2 $fileName


#f1="111014_lane_2_4h-sample-1_1.fastq-se_cufflinks"
#f2="111014_lane_2_4h-sample-1_1.fastq-se_cufflinks_human"

#ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/mix_comb_stats.rb $f1 $f2 $fileName

#f1="111014_lane_2_4h-sample-2_1.fastq-se_cufflinks"
#f2="111014_lane_2_4h-sample-2_1.fastq-se_cufflinks_human"

#ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/mix_comb_stats.rb $f1 $f2 $fileName


#f1="111014_lane_2_4h-sample-3_1.fastq-se_cufflinks"
#f2="111014_lane_2_4h-sample-3_1.fastq-se_cufflinks_human"

#ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/mix_comb_stats.rb $f1 $f2 $fileName

echo "done calculation"
