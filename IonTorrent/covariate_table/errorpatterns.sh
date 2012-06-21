#!/bin/bash
# -cwd

fasta=$1
sam=$2
dir=$3

/nfs/apps/matlab/current/bin/matlab -singleCompThread -nodesktop -nosplash -nodisplay -nojvm -r "cd $dir; errorpatterns('$fasta', '$sam'); quit"
