if [[ $chain != "0" ]]; then ## call realign
    
    OUTDIR=$output.bam_refine

    if [ ! -d $OUTDIR ]; then
	mkdir $OUTDIR
    fi  
    status=$OUTDIR"/realign.status"
    if [ -e $status ] ; then
	rm -f $status
    fi

    touch $status
    mkdir -p $OUTDIR/logs/
    qmem=5 # default
    heapm=4
    for (( i=1; i<=24; i++))
      do 
      if [[ $i -lt 7 ]]; then  # mem=8 for large chr
	  qmem=8
	  heapm=7
      fi
      
      g=`basename $output.bam | sed 's/\//_/g'`
	if [[ $AUTO == "" ]];then
              cmd="qsub -N realign.$i.$g -l mem=${qmem}G,time=55:: -o $OUTDIR/logs/realign.$i.o -e $OUTDIR/logs/realign.$i.e ${BPATH}/gatk_realign_atomic.sh -I $output.bam -o $OUTDIR  -g $setting -L $i -c $status -m $heapm "
  	else
	      cmd="qsub -N realign.$i.$g -l mem=${qmem}G,time=55:: -o $OUTDIR/logs/realign.$i.o -e $OUTDIR/logs/realign.$i.e ${BPATH}/gatk_realign_atomic.sh -I $output.bam -o $OUTDIR  -g $setting -L $i -c $status -m $heapm -A AUTO"
	fi
      echo $cmd
      $cmd
    done


###     $QSUB -l mem=8G,time=58::  $BPATH/gatk_realign_all.sh   -g $setting -I $output.bam 


fi
