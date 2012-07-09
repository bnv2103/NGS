#!/bin/sh
#$ -cwd

prefix=$1

rm homoref.txt het.txt homonull.txt
for i in `seq 1 100`
do
	thresh=`echo "$i/100" | bc -l`
	awk -v d=$thresh '{ if($2==0) { if($5<d&&$4>=$6) trueaa++; else aab++; } else { if($5>=d||$4<$6) falseaa++; else baa++; } } END{aasens = trueaa/(trueaa+aab); aaspc = falseaa/(falseaa+baa); print 1-aaspc, "\t", aasens;}' HMM/output.sta >> homoref.txt
	awk -v d=$thresh '{ if($2==1||$2==2) { if($5>=d) trueab++; else abh++; } else { if($5<d) falseab++; else hab++;} }          END{absens = trueab/(trueab+abh); abspc = falseab/(falseab+hab); print 1-abspc, "\t", absens;}' HMM/output.sta >> het.txt
	awk -v d=$thresh '{ if($2==3) { if($5<d&&$6>=$4) truebb++; else bba++; } else { if($5>=d||$6<$4) falsebb++; else abb++; } } END{bbsens = truebb/(truebb+bba); bbspc = falsebb/(falsebb+abb); print 1-bbspc, "\t", bbsens;}' HMM/output.sta >> homonull.txt
done

# Code to get list of all prefix soft clip lengths
grep Start ${prefix}.sorted.sam | cut -f 6 | cut -d'S' -f 1 | grep -v [A-Z]

# Code to obtain set of input read positions and their haplotype assignments
grep -v ^@ ${prefix}.sorted.sam  | cut -f 1,4 | cut -d ':' -f 3 > ${prefix}.obs
paste ${prefix}.obs HMM/output.obs | awk '{if($2!=$4) printf "reads %d and %d do not tally\n",$2,$4; else if($1!=$3) printf "haplotypes %d and %d for reads %d and %d do not tally\n", $1, $3, $2, $4; else {printf "haplotypes %d and %d for reads %d and %d tally\n", $1, $3, $2, $4; sum++;} } END{printf "Total %d tally from %d: %f\n", sum, NR, (sum/NR)*100;}'

# Code to get ratio of correct genotype assignments
awk '{ if($2==0) {tote++;if($3==0) sume++;} else if($2==1) {totk++;if($3==1) sumk++;} else if($2==2) {totn++;if($3==1) sumn++;} else if($2==3) {toth++;if($3==2) sumh++;} else print "known value: ", $2;} END{print sume,sume/tote;print sumk,sumk/totk;print sumn,sumn/totn;print sumh,sumh/toth;}' HMM/output.sta

