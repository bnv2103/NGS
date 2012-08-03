#!/bin/sh
#$ -cwd
## 10G,1::

id=$1

goby 8g sort -d tmp -s 10000000 --output ${id} ${id}

## --format=goby --goby-output=Sample01 Sample01.fq 

