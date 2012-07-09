#!/bin/sh
#$ -cwd
## 20G,4::

chr=$1

cmetindex -d $chr -k 15 -b 15

