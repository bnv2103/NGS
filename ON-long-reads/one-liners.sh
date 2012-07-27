#!/bin/sh
#$ -cwd
#############################################################################################################
######## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ########
######### NEVER EVER EVER RUN THIS SCRIPT BY ITSELF. PICK UP THE CODE INDIVIDUALLY AND RUN IT ALONE #########
######## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ########
#############################################################################################################

prefix=$1

# Code to extract phased snp info in desired format
grep -v ^# reference/CEU.trio.2010_03.genotypes.vcf | awk '{if(match($12,/[0-1]\|[0-1]:*/)) {split($12,gt,":"); split(gt[1],val,"|"); if(val[1]==0) mom=$4; else mom=$5; if(val[2]==0) dad=$4; else dad=$5; if(match($3,"rs")) known=1; else known=0; print $1, $2, $4, $5, val[1], val[2], mom, dad, known }}' > reference/snp.list

grep -v ^# reference/CEU.trio.2010_07.indel.genotypes.vcf | awk '{split($12,gt,":"); split(gt[1],val,"/"); split($5,alt,","); if(val[1]==0) mom=$4; else mom=alt[val[1]]; if(val[2]==0) dad=$4; else dad=alt[val[2]]; print $1, $2, $4, $5, val[1], val[2], mom, dad }' > reference/indel.list

# Code to find sensitivity to accuracy ratio based on maximum value of genotype - single point in plot
# Maybe split for known and novel hets later?
awk '{if($2==1||$2==2) {tot++;if($3==1) sum++;} else {pot++; if($3!=1) pum++;}} END{sens = sum/tot; spc = pum/pot; print sens, "\t", spc;}' HMM/output.sta > max.txt

#############################################################################################################
######## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ########
######### NEVER EVER EVER RUN THIS SCRIPT BY ITSELF. PICK UP THE CODE INDIVIDUALLY AND RUN IT ALONE #########
######## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ########
#############################################################################################################

# Code to obtain sensitivity to accuracy ratio for different cut-off for heterozygous snps
# Maybe split for known and novel hets later?
rm homoref.txt het.txt homonull.txt
for i in `seq 1 100`
do
	thresh=`echo "$i/100" | bc -l`
	# awk -v d=$thresh '{ if($2==0) { if($5<d&&$4>=$6) trueaa++; else aab++; } else { if($5>=d||$4<$6) falseaa++; else baa++; } } END{aasens = trueaa/(trueaa+aab); aaspc = falseaa/(falseaa+baa); print 1-aaspc, "\t", aasens;}' HMM/output.sta >> homoref.txt
	awk -v d=$thresh '{ if($2==1||$2==2) { if($5>=d) trueab++; else abh++; } else { if($5<d) falseab++; else hab++;} } END{absens = trueab/(trueab+abh); abspc = falseab/(falseab+hab); print 1-abspc, "\t", absens;}' HMM/output.sta >> het.txt
	# awk -v d=$thresh '{ if($2==3) { if($5<d&&$6>=$4) truebb++; else bba++; } else { if($5>=d||$6<$4) falsebb++; else abb++; } } END{bbsens = truebb/(truebb+bba); bbspc = falsebb/(falsebb+abb); print 1-bbspc, "\t", bbsens;}' HMM/output.sta >> homonull.txt
done

#############################################################################################################
######## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ########
######### NEVER EVER EVER RUN THIS SCRIPT BY ITSELF. PICK UP THE CODE INDIVIDUALLY AND RUN IT ALONE #########
######## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ########
#############################################################################################################

# Code to find area under curve of sensitivity to accuracy plot
# calculated using G=2AUC-1 where, G = 1 - sum(X(k)-X(k-1)*Y(k)+Y(k-1)), k= {2,..,n}
awk '{a[NR]=$1;b[NR]=$2;} END{for(i=NR-1;i>=0;i--) {sum += ((a[i]-a[i+1])*(b[i]+b[i+1]))} print 1-sum/2;}' het.txt > auc.txt

# Code to find point of maximum distance on the plot from a (random) diagonal
awk 'BEGIN{dist=0;x=0;y=0;} {if($2-$1>dist) {dist = $2-$1;y=$1;x=$2;}} END{print x,y,dist;}' het.txt > diag.txt

exit

#############################################################################################################
######## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ########
######### NEVER EVER EVER RUN THIS SCRIPT BY ITSELF. PICK UP THE CODE INDIVIDUALLY AND RUN IT ALONE #########
######## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ########
#############################################################################################################

# Code to get list of all prefix soft clip lengths
grep Start ${prefix}.sorted.sam | cut -f 6 | cut -d'S' -f 1 | grep -v [A-Z]

# Code to obtain set of input read positions and their haplotype assignments
grep -v ^@ ${prefix}.sorted.sam  | cut -f 1,4 | cut -d ':' -f 3 > ${prefix}.obs
paste ${prefix}.obs HMM/output.obs | awk '{if($2!=$4) printf "reads %d and %d do not tally\n",$2,$4; else if($1!=$3) printf "haplotypes %d and %d for reads %d and %d do not tally\n", $1, $3, $2, $4; else {printf "haplotypes %d and %d for reads %d and %d tally\n", $1, $3, $2, $4; sum++;} } END{printf "Total %d tally from %d: %f\n", sum, NR, (sum/NR)*100;}'

# Code to get ratio of correct genotype assignments
awk '{ if($2==0) {tote++;if($3==0) sume++;} else if($2==1) {totk++;if($3==1) sumk++;} else if($2==2) {totn++;if($3==1) sumn++;} else if($2==3) {toth++;if($3==2) sumh++;} else print "known value: ", $2;} END{print sume,sume/tote;print sumk,sumk/totk;print sumn,sumn/totn;print sumh,sumh/toth;}' HMM/output.sta

#############################################################################################################
######## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ########
######### NEVER EVER EVER RUN THIS SCRIPT BY ITSELF. PICK UP THE CODE INDIVIDUALLY AND RUN IT ALONE #########
######## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ########
#############################################################################################################

