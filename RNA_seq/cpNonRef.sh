#!/bin/bash
#$ -cwd

 dirName="nonRef_to_release"
 if [ ! -d $dirName ]; then mkdir -p $dirName; fi

 # for f in *cufflinks; do  sh /ifs/scratch/c2b2/ngs_lab/xs2182/code/cpIsoforms.sh  $f $dirName; done
 for f in *cufflinks; do  sh /ifs/scratch/c2b2/ngs_lab/xs2182/code/cpIsoforms_nonRef.sh  $f $dirName; done
 cd $dirName

 
 Rscript /ifs/scratch/c2b2/ngs_lab/xs2182/code/cleanupsampleName.R 

#  rm *nonRef
