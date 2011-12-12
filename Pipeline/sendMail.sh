#!/bin/sh
#$ -S /bin/sh
#$ -cwd


USAGE="Usage: $0 -s Subject  -t TO[,TO]  -m Message_file -h"

#email subject
SUBJECT='EMPTY SUBJECT'

#email TO
# EMAIL="sz2317@c2b2.columbia.edu,xs2182@c2b2.columbia.edu,yshen@c2b2.columbia.edu,oc2121@c2b2.columbia.edu"
EMAIL="sz2317@c2b2.columbia.edu,yshen@c2b2.columbia.edu,oc2121@c2b2.columbia.edu"

#email message
MESSAGE_FILE=""


while getopts t:s:m:h opt
  do
  case "$opt" in
      t) EMAIL="$OPTARG";;
      m) MESSAGE_FILE="$OPTARG";;
      s) SUBJECT="$OPTARG";;
      h) echo $USAGE; exit 1
  esac
done

if [[ $MESSAGE_FILE == "" ]]
    then
    echo $USAGE
    exit 1
fi


/bin/mail -s "$SUBJECT" $EMAIL < $MESSAGE_FILE

