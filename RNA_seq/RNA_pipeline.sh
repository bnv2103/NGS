#!/bin/bash
#$ -cwd

# genome = mouse, human, both
genome=$1
echo $genome
# fastq files, if it is Paired-end, then 2 fastq files are required
f1=$2
f2=$3

if [[ $genome == "mouse" ]]
    then
    setting_genome="/ifs/scratch/c2b2/ngs_lab/xs2182/code/global_setting_mouse.sh"
fi

if [[ $genome == "human" ]]
    then
    setting_genome="/ifs/scratch/c2b2/ngs_lab/xs2182/code/global_setting_human.sh"
fi

script_SE="/ifs/scratch/c2b2/ngs_lab/xs2182/code/pipeline_cufflink-ref.sh"
script_PE="/ifs/scratch/c2b2/ngs_lab/xs2182/code/pipeline_cufflink-ref-PE.sh"

if [ ! -d logs ]; then mkdir -p logs; fi
# print title for statistical summary
ruby /ifs/scratch/c2b2/ngs_lab/xs2182/code/printTitle.rb "summary.csv"

f1_base=`basename $f1`

if [[ $f2 == "" ]]; then
    qsub -o logs/pipe.$f1.o -e logs/pipe.$f1.e -pe smp 4 -R y -l mem=2G,time=100:: -N do.$f1 $script_SE  -i $f1 -s $setting_genome -n 4 -o "cufflinks/"$f1_base"_"$genome"_se_cufflinks"
    else
    qsub -o logs/pipe.$f1.o -e logs/pipe.$f1.e -pe smp 4 -R y -l mem=2G,time=100:: -N do.$f1 $script_PE  -f $f1 -r $f2 -s $setting_genome -n 4 -o "cufflinks/"$f1_base"_"$genome"_pe_cufflinks"
fi


