#!/bin/sh
#$ -cwd
## 10G,8::

id1=$1
id2=$2

goby 8g ftc --quality-encoding Sanger ${id1} ${id2}

