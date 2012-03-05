#!/bin/bash
#$ -S /bin/sh
#$ -cwd

fqdir=$1
runid=$2

SampleSheets="/ifs/data/c2b2/ngs_lab/ngs/status/SampleSheets/"
DIR="/ifs/scratch/c2b2/ngs_lab/ngs/Projects/"
APP=""

RESOURCES="/ifs/data/c2b2/ngs_lab/ngs/resources/"
EXOMEBASE="/ifs/data/c2b2/ngs_lab/ngs/code/NGS/Exome/"
RNABASE="/ifs/data/c2b2/ngs_lab/ngs/code/NGS/RNA_seq/"

for fq in `ls $fqdir/*_1.fastq`
do
	#E.x. fq=/ifs/scratch/c2b2/ngs_lab/ngs/Fastq/120217_SN828_0119_BD07NDACXX/demultiplex/120217_SN828_0119_BD07NDACXX_lane1_2_1.fastq
	sampleid=`basename $fq | cut -f6 -d '_' ` ;
	lane=`basename $fq | cut -f5 -d '_' |sed 's/lane//' ` ; 
	barcode=`basename $fq | cut -f7 -d '_' ` ;

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
continue;

        if [[ $app =~ "rna-seq" ]];
        then
		APP="RNA-seq/"
	elif  [[ $app =~ "exome" ]];
        then
		APP="Exome-seq/"
	else
		continue;
	fi
	
	if [ ! -d $DIR/$APP/$projectid ];then mkdir -p $DIR/$APP/$projectid; fi
	if [ ! -d $DIR/$APP/$projectid/$runid ];then mkdir -p $DIR/$APP/$projectid/$runid; fi
        if [ ! -d $DIR/$APP/$projectid/$runid/fastq ];then mkdir -p $DIR/$APP/$projectid/$runid/fastq; fi

	#make links to  fastq files in destination dir
	ln -s $fq $DIR/$APP/$projectid/$runid/fastq/
	fq_3=`echo $fq | sed 's/_1.fastq/_3.fastq/' `
	if [ -e $fq_3 ];then ln -s $fq_3 $DIR/$APP/$projectid/$runid/fastq/; fi 	#if the _3.fastq exists, link it also
	base_fq=`basename $fq`
	ln_fq="$DIR/$APP/$projectid/$runid/fastq/$base_fq"
	if [ -e $fq_3 ]; then
		ln_fq_3=`echo $ln_fq | sed 's/_1.fastq/_3.fastq/' `
	else
		ln_fq_3=""
	fi

	if [ ! -d $DIR/$APP/$projectid/$runid/logs ];then mkdir -p $DIR/$APP/$projectid/$runid/logs; fi
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
		if [[ $flag == 0 ]];then
			echo "ERROR: Organism Unknown"
		fi

	elif [[ $app =~ "exome" ]];
	then
		#initiate exome-seq pipeline
		if [ ! -d $DIR/$APP/$projectid/$runid/mapping ]; then mkdir -p $DIR/$APP/$projectid/$runid/mapping; fi
		if [ ! -e global_setting_b37.sh ];then
		        cp $EXOMEBASE/global_setting_b37.sh $DIR/$APP/$projectid/$runid/.
			if [[ $capture =~ "mouse" ]]; then
	                        echo -e "export ExonFile="/ifs/data/c2b2/ngs_lab/ngs/resources/Agilent/SureSelect_All_Exon_V1_with_annotation.Mouse.bed.mod"">> $DIR/$APP/$projectid/$runid/global_setting_b37.sh
		        elif [[ $capture =~ "44mb" ]]; then
		                echo -e "export ExonFile="/ifs/data/c2b2/ngs_lab/ngs/resources/Agilent/SureSelect_All_Exon_V2_with_annotation.hg19.bed.mod"" >> $DIR/$APP/$projectid/$runid/global_setting_b37.sh
		        fi
		fi

		#Call Mapping
		sh $EXOMEBASE/do_mapping.sh $ln_fq $DIR/$APP/$projectid/$runid/global_setting_b37.sh $projectid
	fi

done

