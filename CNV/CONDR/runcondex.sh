#!/bin/sh
#$ -S /bin/sh
#$ -cwd
#$ -l h_vmem=8G,time=8::


CONDEX="/ifs/scratch/c2b2/ngs_lab/sz2317/softwares/CONDR/CONDEX.jar"
DIR="/ifs/scratch/c2b2/ngs_lab/sz2317/softwares/CONDR/test/referenceSamples/"
OUT=$DIR"/outcondr/"
TEMP="$OUT/temp"

if [[ ! -e $TEMP ]]; 
then
	mkdir $TEMP
fi

JAVA="java -Xmx500M -Xms500m -Djava.io.tmpdir="${TEMP}



$JAVA -jar $CONDEX -c "10-10" -e $DIR"sample.exon.chr10" -t -b $DIR"sample1.chr10.exon",$DIR"sample10.chr10.exon",$DIR"sample11.chr10.exon",$DIR"sample12.chr10.exon",$DIR"sample13.chr10.exon",$DIR"sample14.chr10.exon",$DIR"sample15.chr10.exon",$DIR"sample16.chr10.exon",$DIR"sample2.chr10.exon",$DIR"sample3.chr10.exon",$DIR"sample4.chr10.exon",$DIR"sample5.chr10.exon",$DIR"sample6.chr10.exon",$DIR"sample7.chr10.exon",$DIR"sample8.chr10.exon",$DIR"sample9.chr10.exon" -p $DIR"SimulationParameterFile.04102011.length200000.rate500000" -o $OUT"sample.chr10.result" -threshold 0 


