#!/bin/sh
#$ -cwd

prefix=$1
inp=$2

# Code to find sensitivity to accuracy ratio based on maximum value of genotype - single point in plot
# Maybe split for known and novel hets later?
echo "Sensitivity and accuracy(1-specificity) based on highest genotype value" > ${prefix}.a
awk '{if($2==1||$2==2) {tot++;if($3==1) sum++;} else {pot++; if($3!=1) pum++;}} END{sens = sum/tot; spc = pum/pot; print 1-sens, "\t", spc;}' ${prefix}.sta >> ${prefix}.a
echo >> ${prefix}.a

# Code to obtain sensitivity to accuracy ratio for different cut-off for heterozygous snps
# Maybe split for known and novel hets later?
rm ${prefix}.b
for i in `seq 1 100`
do
	thresh=`echo "$i/100" | bc -l`
	awk -v d=$thresh '{ if($2==1||$2==2) { if($5>=d) trueab++; else abh++; } else { if($5<d) falseab++; else hab++;} } END{absens = trueab/(trueab+abh); abspc = falseab/(falseab+hab); print 1-abspc, "\t", absens;}' ${prefix}.sta >> ${prefix}.b
done

# Code to find area under curve of sensitivity to accuracy plot
# calculated using G=2AUC-1 where, G = 1 - sum(X(k)-X(k-1)*Y(k)+Y(k-1)), k= {2,..,n}
echo "Area under curve" >> ${prefix}.a
awk '{a[NR]=$1;b[NR]=$2;} END{for(i=NR-1;i>=0;i--) {sum += ((a[i]-a[i+1])*(b[i]+b[i+1]))} print 1-sum/2;}' ${prefix}.b >> ${prefix}.a
echo >> ${prefix}.a

# Code to find point of maximum distance on the plot from a (random) diagonal
echo "Farthest point from diagonal" >> ${prefix}.a
awk 'BEGIN{dist=0;x=0;y=0;} {if($2-$1>dist) {dist = $2-$1;y=$1;x=$2;}} END{print x,y,dist;}' ${prefix}.b >> ${prefix}.a
echo >> ${prefix}.a

# Code to obtain set of input read positions and their haplotype assignments
paste ${inp}.ibs ${prefix}.obs | awk '{if($2!=$4) printf "reads %d and %d do not tally\n",$2,$4; else if($1!=$3) printf "haplotypes %d and %d for reads %d and %d do not tally\n", $1, $3, $2, $4; else {printf "haplotypes %d and %d for reads %d and %d tally\n", $1, $3, $2, $4; sum++;} } END{printf "Total %d tally from %d: %f\n", sum, NR, (sum/NR)*100;}' >> ${prefix}.a

# Maybe I'll combine these in the super master pipeline script to plot for all simulation parameters in single plot
cp base.gnu ${prefix}.gnu
echo "set output '${prefix}.png'" >> ${prefix}.gnu
echo "plot '${prefix}.b' w l ls 1" >> ${prefix}.gnu
gnuplot ${prefix}.gnu

