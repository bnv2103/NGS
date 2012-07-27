#!/bin/bash
#$ -S /bin/sh
#$ -cwd

fqdir=$1
runid=$2
setting=$3

if [[ $fqdir == "" || $runid == ""  || $setting == "" ]]
    then
    echo "Error: Missing Argument :: $0 $1 $2 $3"
    echo "USAGE: $0 path/fqdir runID path/pipeline_global_settings_file "
    exit 1
fi

. $setting

DIR="/ifs/scratch/c2b2/ngs_lab/ngs/Projects/"
APP=""

for fq in `ls $fqdir/*_1.fastq`
do
	#E.x. fq=/ifs/scratch/c2b2/ngs_lab/ngs/Fastq/120217_SN828_0119_BD07NDACXX/demultiplex/120217_SN828_0119_BD07NDACXX_lane1_2_1.fastq
	sampleid=`basename $fq | cut -f6 -d '_' ` ;
	lane=`basename $fq | cut -f5 -d '_' |sed 's/lane//' ` ; 
	barcode=`basename $fq | cut -f7 -d '_' ` ;
#	rundate=`echo $runid | cut -f1,3 -d'_' ` 

	#search for sample in tsv file, return organism, application, capture, reads, service, projectID
	#if any fields have blanks they are replaced with "_" for ease. Fields expected to have blanks are organism, capture, reads, application

	arr=(`awk -F '\t' '{sample=$5;  gsub(/ /, "-",sample); gsub(/_/, "-",sample); gsub(/\//, "-",sample);   if ( $2 == '${lane}' &&  sample == "'${sampleid}'" && $7 == "'${barcode}'"  ) {org =$6; app= $8; cap = $9; reads=$10; gsub(/ /, "_",org); gsub(/ /, "_",app); gsub(/ /, "_",cap); gsub(/ /, "_",reads);  print tolower(org) ,tolower(app) , tolower(cap), reads, tolower($12), $13; }}' $SampleSheets/$runid.tsv |  tr " " "\n" `)


	if [[ ${#arr[@]} == 0 ]];	
	then
		continue;	#if no line retrieved then check next file
	fi

	#check if we need to do analytics or not
	service=`echo ${arr[4]}| tr '[:upper:]' '[:lower:]' `
	if [[ $service =~ "self" ]];
	then
		continue; 	#if service=self then we donot need to analyse
	fi

	organism=${arr[0]} 	 	#(Human, Mouse, Xenograft etc.)
	app=${arr[1]}			#(Whole exome, RNA-seq, DNA, CHIP-seq, exome)
	capture=${arr[2]}		#(Agilent 44Mb, N/A, Agilent 50Mb)
	reads=${arr[3]}			#(SE,PE, SE100 etc)
	projectid=${arr[5]}

	echo "$lane $sampleid $barcode $projectid"


        if [[ $app =~ "rna-seq" ]];
        then
		APP="RNA-seq/"
	elif  [[ $app =~ "exome" ]];
        then
		APP="Exome-seq/"
	else
		## Fire QC on the Fastq
	        if [ ! -d $fqdir/QC ];then mkdir -p $fqdir/QC ; fi
	        if [ ! -d $fqdir/logs ];then mkdir -p $fqdir/logs ; fi
		j=`basename $fq`
		fq3=`echo $fq | sed 's/_1.fastq/_3.fastq' `
		scriptfile="$fqdir/QC/$j.runQC.sh"
		echo '#!/bin/bash ' > $scriptfile 
		echo '#$ -cwd ' >> $scriptfile
		echo " ruby $UTILS/fastq_QCplot.rb $fqdir/QC/$j.QC $fq $fq3  & " >> $scriptfile
	        echo " sh  $UTILS/picard_LibraryComplexity.sh -i $fq -p  $fq3 -o $fqdir/QC/$j.LibComplexity -m 7 -s $j " >> $scriptfile
		echo " wait " >>  $scriptfile
		## Finally GZip the demuxed file.
		echo " qsub -o $fqdir/logs/zip.$j.o -e $fqdir/logs/zip.$j.e -N Zip.$j -l mem=512M,time=2::  $NGSSHELL/do_gzip.sh $fq $fq3 " >> $scriptfile
		CMD="qsub -o $fqdir/logs/QC.$j.o -e $fqdir/logs/QC.$j.e -N QC.$j -l mem=10G,time=12:: $scriptfile "
		echo $CMD
		$CMD
		continue;
	fi
	
	if [ ! -d $DIR/$APP/$projectid ];then mkdir -p $DIR/$APP/$projectid; fi
	if [ ! -d $DIR/$APP/$projectid/$runid ];then mkdir -p $DIR/$APP/$projectid/$runid; fi
        if [ ! -d $DIR/$APP/$projectid/$runid/fastq ];then mkdir -p $DIR/$APP/$projectid/$runid/fastq; fi
	if [ ! -d $DIR/$APP/$projectid/$runid/logs ];then mkdir -p $DIR/$APP/$projectid/$runid/logs; fi

	#move fastq files to destination dir
	mv  $fq $DIR/$APP/$projectid/$runid/fastq/
	fq_3=`echo $fq | sed 's/_1.fastq/_3.fastq/' `
	if [ -e $fq_3 ];then mv $fq_3 $DIR/$APP/$projectid/$runid/fastq/; fi 	#if the _3.fastq exists, link it also
	base_fq=`basename $fq`
	ln_fq="$DIR/$APP/$projectid/$runid/fastq/$base_fq"
	ln_fq_3=`echo $ln_fq | sed 's/_1.fastq/_3.fastq/' `

	## Fire QC on the Fastq for RNA-seq & Exome-Seq 
        if [ ! -d $fqdir/QC ];then mkdir -p $fqdir/QC ; fi
        if [ ! -d $fqdir/logs ];then mkdir -p $fqdir/logs ; fi
	scriptfile="$fqdir/QC/$base_fq.runQC.sh"
		echo '#!/bin/bash ' > $scriptfile 
		echo '#$ -cwd ' >> $scriptfile
		echo " ruby $UTILS/fastq_QCplot.rb $fqdir/QC/$base_fq.QC $ln_fq $ln_fq_3  & " >> $scriptfile
	        echo " sh  $UTILS/picard_LibraryComplexity.sh -i $ln_fq  -p $ln_fq_3 -o $fqdir/QC/$base_fq.LibComplexity -m 7 -s $base_fq " >> $scriptfile
		echo " wait " >>  $scriptfile
		## Finally GZip the demuxed file.
		CMD="qsub -o $fqdir/logs/QC.$base_fq.o -e $fqdir/logs/QC.$base_fq.e -N QC.$base_fq -l mem=10G,time=12:: $scriptfile "
		echo $CMD
		$CMD
#        QCDir="$fqdir/QC/$lane_$sampleid/"
#	if [ ! -d $QCDir ];then mkdir -p $QCDir ; fi
	if [ ! -e $ln_fq_3 ]; then
		ln_fq_3=""
	fi
	cd $DIR/$APP/$projectid/$runid

	#Initiate Pipelines according to app!
	if [[ $app =~ "rna-seq" ]];
	then
		#initiate RNA-seq pipeline ;send information organism and se & pe files.
		flag=0
		if [[ $organism =~ "human" ]];then
		     sh $RNABASE/RNA_pipeline.sh human $ln_fq $ln_fq_3
			flag=1
		fi
		if  [[ $organism =~ "mouse" ]];then
        	     sh $RNABASE/RNA_pipeline.sh mouse $ln_fq $ln_fq_3
			flag=1
		fi
		if  [[ $organism =~ "xenograft" ]];then
        	     sh $RNABASE/RNA_pipeline.sh human $ln_fq $ln_fq_3
	             sh $RNABASE/RNA_pipeline.sh mouse $ln_fq $ln_fq_3
			flag=1
		fi
                if  [[ $organism =~ "rat" ]];then
                     sh $RNABASE/RNA_pipeline.sh rat $ln_fq $ln_fq_3
                        flag=1
                fi
                if  [[ $organism =~ "fruitfly" ]];then
                     sh $RNABASE/RNA_pipeline.sh fruitfly $ln_fq $ln_fq_3
                        flag=1
                fi
#                if  [[ $organism =~ "yeast" ]];then
#                     sh $RNABASE/RNA_pipeline.sh yeast $ln_fq $ln_fq_3
#                        flag=1
#                fi

		if [[ $flag == 0 ]];then
			echo "ERROR: Organism Unknown"
		fi

	elif [[ $app =~ "exome" ]];
	then
		#initiate exome-seq pipeline
		EXOMEBASE_1=$EXOMEBASE
		if [ ! -d $DIR/$APP/$projectid/$runid/mapping ]; then mkdir -p $DIR/$APP/$projectid/$runid/mapping; fi
		if [ ! -e global_setting_b37.sh ];then
		        cp $EXOMEBASE_1/global_setting_b37.sh $DIR/$APP/$projectid/$runid/global_setting.sh
			if [[ $capture =~ "mouse" ]]; then
				EXOMEBASE_1="/ifs/data/c2b2/ngs_lab/ngs/code/NGS/Exome_Mouse/";
				cp $EXOMEBASE_1/global_setting.sh $DIR/$APP/$projectid/$runid/global_setting.sh
		        elif [[ $capture =~ "44mb" ]]; then
		                echo -e "export ExonFile="/ifs/data/c2b2/ngs_lab/ngs/resources/Agilent/SureSelect_All_Exon_V2_with_annotation.hg19.bed.mod"" >> $DIR/$APP/$projectid/$runid/global_setting.sh
                        elif [[ $capture =~ "51mb" ]]; then
                                echo -e "export ExonFile="/ifs/data/c2b2/ngs_lab/ngs/resources/Agilent/SureSelect_All_Exon_V4_hg19.bed.sorted.bed"" >> $DIR/$APP/$projectid/$runid/global_setting.sh
		        fi
		fi
		#Call Mapping
		CMD="sh $EXOMEBASE_1/do_mapping.sh $ln_fq $DIR/$APP/$projectid/$runid/global_setting.sh $projectid "
		echo $CMD
		$CMD
	fi
done

