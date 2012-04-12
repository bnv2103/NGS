#!/bin/bash
#$ -cwd

dir=$1
filesDir=$2

cd /ifs/data/shares/deepsequencing/solid/public_html/

mkdir $dir
chmod 775 $dir
cd $dir

for f in $filesDir"/*"; do
cp $f "/ifs/data/shares/deepsequencing/solid/public_html/"$dir
done


cd ..
sh /ifs/data/c2b2/ngs_lab/ngs/code/NGS/Pipeline/make_release_index-html.sh $dir