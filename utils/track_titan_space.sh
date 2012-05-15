#!/bin/bash
#$ -cwd

suff=` date +%N`
dir=`pwd`

to_email="oc2121@c2b2.columbia.edu,sz2317@c2b2.columbia.edu,xs2182@c2b2.columbia.edu,yshen@c2b2.columbia.edu,wsd2102@c2b2.columbia.edu" 

total=`df -h /ifs/scratch/c2b2/ngs_lab/ngs/ |tail -1 | sed 's/ \+ /\t/g' |cut -f2`
used=`df -h /ifs/scratch/c2b2/ngs_lab/ngs/ |tail -1 | sed 's/ \+ /\t/g'  |cut -f3`
avail=`df -h /ifs/scratch/c2b2/ngs_lab/ngs/ |tail -1 | sed 's/ \+ /\t/g' |cut -f4`

if [[ $avail =~ "T" ]];then
	tera=`echo $avail | sed 's/T//'`
	avail_byte=`echo "$tera*1024*1024*1024*1024.0" | bc  -l | awk -F '.' '{ print $1; exit; }'`
elif  [[ $avail =~ "G" ]];then
        tera=`echo $avail | sed 's/T//'`
        avail_byte=`echo "$tera*1024*1024*1024.0" | bc  -l | awk -F '.' '{ print $1; exit; }'`
elif  [[ $avail =~ "M" ]];then
        tera=`echo $avail | sed 's/T//'`
        avail_byte=`echo "$tera*1024*1024.0" | bc  -l | awk -F '.' '{ print $1; exit; }'`
fi

threshold=` echo "4*1024*1024*1024*1024.0" | bc  -l | awk -F '.' '{ print $1; exit; }'`

if [[ $avail_byte -lt $threshold ]];then
	echo "NGS Disk Space Alarm." > $dir/mailbody_$suff.txt	
	echo "Available disk space : $avail " >> $dir/mailbody_$suff.txt
        echo "Used disk space : $used " >> $dir/mailbody_$suff.txt
        echo "Total disk space : $total " >> $dir/mailbody_$suff.txt
	sh /ifs/data/c2b2/ngs_lab/ngs/code/NGS/Pipeline/sendMail.sh -s "NGS Disk Space Alarm" -t $to_email  -m $dir/mailbody_$suff.txt
	rm $dir/mailbody_$suff.txt
fi

