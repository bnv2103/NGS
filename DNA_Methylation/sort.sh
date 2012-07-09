#!/bin/sh
#$ -cwd
## 2G,1::

id=$1

goby 1g sort --output ${id} ${id}

## --format=goby --goby-output=Sample01 Sample01.fq 

