#!/bin/sh
#$ -cwd

prefix=$1

# Code to get list of all prefix soft clip lengths
grep Start ${prefix}.sorted.sam | cut -f 6 | cut -d'S' -f 1 | grep -v [A-Z]

# Code to obtain set of input read positions and their haplotype assignments
grep -v ^@ ${prefix}.sorted.sam  | cut -f 1,4 | cut -d ':' -f 3 > ${prefix}.obs
paste ${prefix}.obs HMM/output.obs | awk '{if($2!=$4) printf "reads %d and %d do not tally\n",$2,$4; else if($1!=$3) printf "haplotypes %d and %d for reads %d and %d do not tally\n", $1, $3, $2, $4; else {printf "haplotypes %d and %d for reads %d and %d tally\n", $1, $3, $2, $4; sum++;} } END{printf "Total %d tally from %d: %f\n", sum, NR, (sum/NR)*100;}'

