#!/bin/bash
#$ -cwd

runname=$1
outdir=$2

cp $runname/Data/reports/Summary/read*.xml $outdir/
