#!/bin/bash
#$ -cwd


rdir=$1

cd $rdir

flist=`ls`

touch index.html
echo -e "<html><body> <p><ul>" >> index.html


for f in $flist;  do
    s=`du -h $f | head -1 | cut -f1`
    echo -e "<li><a href=$f>$f</a> &nbsp; &nbsp; ($s)</li>" >> index.html
done

echo -e "</ul></p></body></html>" >> index.html

`chmod 664 *`

