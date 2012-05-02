#!/usr/bin/perl

#$ -cwd

## Use as  perl /ifs/scratch/c2b2/ngs_lab/sz2317/runs/ALI_GHARAVI/convert_vcf_exomeannotation.pl infilename [indel] 
# give the indel option for indel  input file  

 use File::Basename;

$vcfin = $ARGV[0];
$outfile = basename($vcfin);

$indel ="";
if ($ARGV[1]=~"indel" || $ARGV[1]=~"1")
{ $indel = "1";}


print "\nOpen _VCF file";

open(vcf, "<" . $vcfin) || die("Could not open _vcf file!");
open(out, ">" . $outfile . ".xls" )|| die("Could not open outfile.xls  file!");

print  "\nBeginning convrsion";

### header of vcf file has 10 fields as follows
# #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  136-04
# 1       871146  .       C       T       923.51  PASS    AB=0.492;AC=1;AF=0.042;AN=24;BaseQRankSum=1.084;DP=1029;Dels=0.00;FS=1.545;HRun=2;HaplotypeScore=5.5593;InbreedingCoeff=-0.0435;MQ=57.01;MQ0=0;MQRankSum=-0.674;QD=15.65;ReadPosRankSum=1.003;SB=-332.61;refseq.chr=1;refseq.codingCoordStr=c.306-6;refseq.end=871146;refseq.haplotypeAlternate=T;refseq.haplotypeReference=C;refseq.inCodingRegion=false;refseq.name=NM_152486;refseq.name2=SAMD11;refseq.positionType=intron;refseq.spliceDist=-6;refseq.spliceInfo=splice-acceptor_-6;refseq.start=871146;refseq.transcriptStrand=+  GT:AD:DP:GQ:PL   0/0:117,0:117:99:0,343,4408

# header for output file is
print out "Chr\tPosition\tRS no\tBase change\tPresent 1000 genome\?\tPresent other DB\?\tPresent EVS\?\tGene Symbol\t";
if ($indel == "1") { print out "Functional Class\t"; }
else{ print out "A\.A\.change\t"};
print out "Splice Site change\tQuality Score\tFilter\tDepth";

my $sampleCount =1;
my @samplelist=();
while(<vcf>) {
	if ($_ =~ m/^#/)
	{
		if ($_ =~ m/^#CHROM/){
			chomp;
			my @header = split("\t");
			$sampleCount = scalar(@header);
			for (my $i=9; $i< $sampleCount; $i++)
			{
				print out "\tGT:$header[$i]\tHet:$header[$i]\tGQ:$header[$i]\tAD:$header[$i]";
			}
			print out "\n";
			next;
		}else{next;}
	}
	else {
	       	chomp;
		my @field = split("\t");
		my @info = split(/;/,$field[7]);	
		my @format = split(":",$field[8]);

		my $buffer="chr";
		$buffer.=$field[0]."\t" ;		#Chr
                $buffer.=$field[1]."\t" ;		#Pos
                $buffer.=  ($field[2]  =~ m/^\./)? "-" : $field[2];	#Rs No
                $buffer.= "\t";
                $buffer.=  $field[3].">". $field[4]."\t";		#Base change
                $buffer.=  ($info[0] =~ m/1KG/ ) ? "YES" : "-";		#Present 1000 genome
                $buffer.= "\t";
                $buffer.=  ($field[2]  =~ m/^\./)? "-" : "dbSNP132";	#Presetn other DB
                $buffer.= "\t";
                $buffer.=  ($field[7]  =~ m/EVS/)? "YES" : "-";    #Present in EVS?
                $buffer.= "\t";

		$i=0;
		if ($field[7] !~ /refseq\.name2/ ) { 
        	       $buffer.= "-\t-\t-\t";
		}
		else{
			while ( ($info[$i] !~ m/^refseq\.name2/ ) &&( $i< (scalar(@info) ) )) {	#Gene symbol
                	        $i++;
                	}
			my $name= (split("=",$info[$i]))[1] ;
			if  ($field[7]  =~ m/^refseq\.name2_/)	# means there are multiple transcripts, then extract all gene symbols.
			{
				$i++;
	        	        while (($info[$i] =~ m/^refseq\.name2_/)  &&( $i< (scalar(@info) ) )) {
					my $gene_sym= (split("=",$info[$i]))[1];
					if ($name !~ m/$gene_sym/)
						{$name.= ",".$gene_sym; }
        	                	$i++;
	                	}
			}
			$buffer.= $name."\t";			#gene Symbol

		    if ($indel == "1") {
			if ($field[7] =~ /refseq\.functionalClass/ ) {   
				$i=0;
	                        while (($info[$i] !~ m/^refseq\.functionalClass/) &&( $i< (scalar(@info) ) )) {
        	                        $i++;
	                        }
                        my $fclass= (split("=",$info[$i]))[1];
                        $name= $fclass;
                        if  ($info[$i] =~ m/^refseq\.functionalClass_/) 
                        {
                                $i++;
                                while (($info[$i] =~ m/^refseq\.functionalClass_/) &&( $i< (scalar(@info) ) )) {
                                        $fclass= (split("=",$info[$i]))[1];
                                        if ($name !~ m/$fclass/ ) {
                                           $name.= ",".$fclass;}
                                        $i++;
                        }}
                        $buffer.= $name."\t";
			}else {  $buffer.="-\t";}
		    }
		    else{
			if ($field[7] =~ /refseq\.changesAA=true/ ) {	# A.A. change
				$i=0;
	                	while( ($info[$i] !~ m/^refseq\.proteinCoordSt/)  &&( $i< scalar(@info) )) {
        	                	$i++;
	                	}
				if ($i >= scalar(@info)){
        		                $buffer.= "-\t";}
                		else{
		                 $buffer.= (split("=",$info[$i]))[1];
				 $buffer.="\t";
				}
			}else {  $buffer.="-\t";}
		     }


	            if( $field[7] =~ /refseq\.spliceDist/ ) {	# Splice site change
			$i=0;
                  	while (($info[$i] !~ m/^refseq\.spliceDist/) &&( $i< (scalar(@info) ) )) {
                        	$i++;
                	}
                	my $splice= (split("=",$info[$i]))[1];
                	$name= $splice;
                	if  ($info[$i] =~ m/^refseq\.spliceDist_/)   #means there are multiple transcripts, then extract all gene symbols.
			{
		                $i++;
        		        while (($info[$i] =~ m/^refseq\.spliceDist_/) &&( $i< (scalar(@info) ) )) {
                		        $splice= (split("=",$info[$i]))[1];
	                	        if ($name !~ m/$splice/ ) {
					   $name.= ",".$splice;}
		                        $i++;
				}
        	       	}
			$buffer.= $name;
                     }else {  $buffer.="-";}
		$buffer.= "\t";
		}
		
		$buffer.= $field[5]."\t";		#QualityScore
		$buffer.= $field[6]."\t";        	#Filter

		$i=0;
		while ($info[$i] !~ m/^DP=/) {				#DP coverage
			$i++;
		}
		$buffer.= ($i >= scalar(@info)) ? "-" : (split("=",$info[$i]))[1] ;

		##Samples genotype
	   my $iter=9;
	   while ($iter < $sampleCount){
		$buffer.= "\t";
		my @data = split(":",$field[$iter]);
                $buffer.=  ($data[0] =~ m/^0\// ) ? $field[3] :  $field[4];	#GT (sample Genotype)
                $buffer.= "/";
                $buffer.=  ($data[0] =~ m/\/0$/ ) ? $field[3] :  $field[4];
                $buffer.= "\t";
                
		if ($data[0] =~ m/0\/1/ ) {				#HET. HOM ?
			$buffer.="HET";}
		elsif ($data[0] =~ m/0\/0/ ) {
                        $buffer.="HOM REF";}
		elsif ($data[0] =~ m/1\/1/ ) {
                        $buffer.="HOM ALT";}
		else{
			$buffer.="-";}

		$buffer.= "\t";
		$i=0; 
                while ($format[$i] !~ m/^GQ/) {                           #GQ - sample genotype quality
                        $i++; }
                $buffer.= ($i >= scalar(@format)) ? "-" : $data[$i];
		$buffer.= "\t";

		$i=0;
		while ($format[$i] !~ m/^AD/) {				#AD Allelic depth for samples
			$i++; }
                $buffer.= ($i >= scalar(@format)) ? "-" : $data[$i];

		$iter++;
	    }
	$buffer.= "\n";
	print out $buffer;
	}
}
close (vcf);
close (out);
print  "\nConversion to exome annotation done\n";

exit;
