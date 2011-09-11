#!/bin/bash
#$ -cwd

config=$1

if [[ $config == "" ]]; then
	config="/ifs/data/c2b2/ngs_lab/ngs/status/Config/ionTorrent.config.sh"
fi

## get files from ion torrent server:
# username and password are read from a config file

ionserver="c2b2ocpd1.c2b2.columbia.edu"

. $config

info=`curl --user $username:$password  --header "Content-Type: application/json"  --location http://$ionserver/rundb/api/v1/experiment`

echo $info
