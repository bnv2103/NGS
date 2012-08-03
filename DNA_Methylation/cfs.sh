#!/bin/sh
#$ -cwd
## 10G,8::

id1=$1
id2=$2

goby 3g cfs `cat groupA.list` -o groupA.stats

goby 3g cfs `cat groupB.list` -o groupB.stats

