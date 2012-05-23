#!/bin/bash
#$ -cwd

### workflow:
# two parameters: setting and outdir

setting=$1
outdir=$2

. $setting
 bam=$outdir"/accepted_hits.bam"

## do cufflinks
cuffout=$bam"_cufflinks_ref"
cuffout2=$bam"_cufflinks_ref-guide"

cmd="cufflinks -o $cuffout --GTF $GENES $bam"
# cmd="cufflinks -o $cuffout --compatible-hits-norm --GTF  $GENES $bam"
echo -e "do cufflinks with ref genes: \n $cmd"
$cmd

cmd2="cufflinks -o $cuffout2 --GTF-guide  $GENES $bam"
echo -e "do cufflinks with ref genes -guide: \n $cmd2"
$cmd2

## reference-guided assembly

# cmd2="cufflinks -o $cuffout2 --GTF-guide  $GENES $bam"
# echo -e "do cufflinks with ref genes -guide: \n $cmd2"
# $cmd2





# if [[ $GENO == "mouse" ]];
#     then
#    qsub -l mem=2G,time=10:: /ifs/scratch/c2b2/ngs_lab/xs2182/code/getBed_mouse.sh $bam $outdir
#     qsub -l mem=2G,time=5:: /ifs/scratch/c2b2/ngs_lab/xs2182/code/getSNPs.sh $outdir "mouse"
# fi
# if [[ $GENO == "human" ]];
#     then
#     qsub -l mem=2G,time=10:: /ifs/scratch/c2b2/ngs_lab/xs2182/code/getBed.sh $bam $outdir
#     qsub -l mem=2G,time=5:: /ifs/scratch/c2b2/ngs_lab/xs2182/code/getSNPs.sh $outdir "human"
# fi
# if [[ $GENO == "rat" ]];
#     then
 #    qsub -l mem=2G,time=10:: /ifs/scratch/c2b2/ngs_lab/xs2182/code/getBed_mouse.sh $bam $outdir
#    qsub -l mem=2G,time=5:: /ifs/scratch/c2b2/ngs_lab/xs2182/code/getSNPs.sh $outdir "rat"
#fi






