

# Code to extract methylation and depth of coverage info from methylation vcf file
grep -v ^# variants/meth.vcf | cut -f 8 | awk -F";" 'BEGIN{printf "LOD\tLOD_SE\tLOD_Z\tFisher\tDelta_MR\tMRA\tMRB\tCGroupA\tCmGroupA\tCGroupB\tCmGroupB\tDP\n";} {for(i=4;i<=15;i++) {split($i,vals,"=");printf vals[2]"\t"} printf "\n";}' > MergedSets/anal.txt

cut -f 12 MergedSets/anal.txt | grep -v DP > MergedSets/depth.txt
cut -f 6,7,5 MergedSets/anal.txt | grep -v MR | grep -v NaN > MergedSets/meth.txt

# Code to get a textual distribution of values of the comparitive methylation analysis of two groups
echo;echo;echo;echo;echo;for i in `seq 1 12`; do cut -f $i MergedSets/anal.txt | awk '{if($0=="NaN") nansum++; else if($0=="-Infinity") infsum++; else if($0<0) {negsum++; neg+= $0} else if($0==0) zerosum++; else if($0>0) {possum++; pos += $0;}} END{if(negsum>0&&possum>0) printf "Nan=%d, Inf=%d, negs=%d, average=%f, pos= %d, average=%f, zero=%d\n", nansum, infsum, negsum, neg/negsum, possum, pos/possum, zerosum; else if(negsum>0) printf "Nan=%d, Inf=%d, negs=%d, average=%f, pos= %d, average=%f, zero=%d\n", nansum, infsum, negsum, neg/negsum, possum, pos, zerosum; else if(possum>0) printf "Nan=%d, Inf=%d, negs=%d, average=%f, pos= %d, average=%f, zero=%d\n", nansum, infsum, negsum, neg, possum, pos/possum, zerosum; else printf "Nan=%d, Inf=%d, negs=%d, average=%f, pos= %d, average=%f, zero=%d\n", nansum, infsum, negsum, neg, possum, pos, zerosum;}' ; done

