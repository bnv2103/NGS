#!/bin/bash
#$ -cwd

# genome = mouse, human, both
genome=$1
# echo $genome
# fastq files, if it is Paired-end, then 2 fastq files are required
f1=$2
f2=$3

#RNABASE="/ifs/data/c2b2/ngs_lab/ngs/code/NGS/RNA_seq/"
RNABASE="/ifs/scratch/c2b2/ngs_lab/xs2182/code/"
RUBY18="/ifs/data/c2b2/ngs_lab/ngs/usr/local/bin/ruby"

if [[ $genome == "mouse" ]];
    then
   setting_genome="$RNABASE/global_setting_mouse.sh"
fi

if [[ $genome == "human" ]];
    then
    setting_genome="$RNABASE/global_setting_human.sh"
fi

if [[ $genome == "rat" ]];
    then
    setting_genome="$RNABASE/global_setting_rat.sh"
fi

if [[ $genome == "fruitfly" ]];
    then
    setting_genome="$RNABASE/global_setting_fruitfly.sh"
fi

if [[ $genome == "yeast" ]];
    then
    setting_genome="$RNABASE/global_setting_yeast.sh"
fi


script_SE="$RNABASE/pipeline_cufflink-ref.sh"
# script_SE="$RNABASE/cufflinkOnly.sh"
script_PE="$RNABASE/pipeline_cufflink-ref-PE.sh"

if [ ! -d  logs ]; then mkdir -p logs; fi
if [ ! -d  cufflinks ]; then mkdir -p cufflinks; fi
if [ ! -d  spikeIn ]; then mkdir -p spikeIn; fi
# print title for statistical summary
# ruby $RNABASE/printTitle.rb "summary.csv"

f1_base=`basename $f1`

if [[ $f2 == "" ]]; then
    ruby $RNABASE/printTitle.rb "summary.csv"
    qsub -o logs/pipe.$f1_base.o -e logs/pipe.$f1_base.e -pe smp 4 -R y -l mem=2G,time=100:: -N do.$f1_base $script_SE  -i $f1 -s $setting_genome -n 4 -o "cufflinks/"$f1_base"_"$genome"_se_cufflinks"
else
    ruby $RNABASE/printTitle_PE.rb "summary.csv"
    qsub -o logs/pipe.$f1_base.o -e logs/pipe.$f1_base.e -pe smp 4 -R y -l mem=4G,time=100:: -N do.$f1_base $script_PE  -f $f1 -r $f2 -s $setting_genome -n 4 -o "cufflinks/"$f1_base"_"$genome"_pe_cufflinks"
fi



qsub -o logs/pipe.$f1.o -e logs/pipe.$f1.e  -l mem=8G,time=20:: -N do.$f1 /ifs/scratch/c2b2/ngs_lab/xs2182/code/pipeline_spikein.sh  -i $f1 -s /ifs/scratch/c2b2/ngs_lab/xs2182/code/global_setting_spike.sh -n 1 -o "spikeIn/"$f1"_spikein";




