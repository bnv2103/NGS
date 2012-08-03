#!/bin/sh
#$ -cwd

qsub -l mem=32G,time=8:: ./methylation_pipeline.sh file.list groupA.list groupB.list mm10

