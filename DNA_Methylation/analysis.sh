#!/bin/sh
#$ -cwd
### mem=6G,time=16

regions=$1
base=$2
range=$3
init=$4
fileA=$5
fileB=$6

ngroups=2

grep -v ^# ../variants/meth.vcf | cut -f 8 | awk -F";" 'BEGIN{printf "LOD\tLOD_SE\tLOD_Z\tFisher\tDelta_MR\tMRA\tMRB\tCGroupA\tCmGroupA\tCGroupB\tCmGroupB\tDP\n";} {for(i=4;i<=15;i++) {split($i,vals,"=");printf vals[2]"\t"} printf "\n";}' > anal.txt

cut -f 12 anal.txt | grep -v DP > depth.txt
cut -f 6,7,5 anal.txt | grep -v MR | grep -v NaN | awk '{print log($1+1)/log(2),"\t",log($2+1)/log(2),"\t",$3}' > meth.txt

if [[ $init -eq 1 ]]
then
	goby 1g annotations-to-counts ${regions}.tsv -o ${regions}
fi

qsub -sync y -l mem=2G,time=2:: -t 1:$range ./counts.sh $base

for i in `seq 1 $range`
do
	goby 1g coverage ${base}${i} -a ${regions} -o ${base}${i}.tsv
done

goby 5g alignment-to-annotation-counts ${base}*.entries --annotation mm10-exons-ensembl-annotation.txt --include-annotation-types gene --compare A/B --groups A=`cat ${fileA} | tr '\n' ','`/B=`cat ${fileB} | tr '\n' ','` --stats stats.tsv

for i in `seq 1 $range`
do
	cut -f 2,10,11,13 ${base}${i}.ann-counts.tsv | grep ENS | awk '{print $1,"\t",log($2+1)/log(2),"\t",log($3+1)/log(2),"\t",$4}' > stat-${base}${i}.txt
done

let "rpkm1=2+(3*$range)+$ngroups+1"
let "rpkm2=$rpkm1+1"
let "fold=$rpkm2+$ngroups+1"
let "foldmag=$fold+1"
let "foldlog=$foldmag+1"

cut -f 1,${rpkm1},${rpkm2},${fold},${foldmag},${foldlog} stats.tsv | grep ^ENS > stat.txt

gnuplot depth.gnu
gnuplot meth.gnu
gnuplot delta-meth.gnu

gnuplot fold-change.gnu
gnuplot fold-change-log.gnu
gnuplot fold-change-mag.gnu
gnuplot RPKM-log.gnu

gnuplot stat-in-count.gnu
gnuplot stat-over-count.gnu
gnuplot stat-rpkm-log.gnu

