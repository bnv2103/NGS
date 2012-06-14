#!/usr/bin/perl

#For Mouse : Convert gene annotation from Annovar back to VCF 
my $fin1 = $ARGV[0];   # variant_function file
my $fin2 = $ARGV[1];   # exonic.variant_function file
my $fin3 = $ARGV[2];   # original VF file, to grab header 

#VCF files assuming chr 1 to 22 and X and Y only
#headerssss
#1       3421959 other fields..
#1       3688030 other fields..

# head -2 mapping/1.bam_refine/VarCalling/list.vcf-files.txt.vcf.dbsnp.annovar.variant_function
# exonic  Rrs1    chr1    9536234 9536234 T       G       chr1    9536234 rs31477593      T       G       26.78   LowQual AB=0.730;AC=1;AF=0.50;AN=2;BaseQRankSum=-1.453;DP=37;Dels=0.00;FS=9.688;HRun=2;HaplotypeScore=0.9967;MQ=52.19;MQ0=0;MQRankSum=-4.497;QD=0.72;ReadPosRankSum=-0.804;SB=-0.01;snp128   GT:AD:DP:GQ:PL  0/1:27,10:37:56.77:57,0,768
# exonic  Rrs1    chr1    9536252 9536252 A       G       chr1    9536252 rs13467800      A       G       105.92  PASS    AB=0.675;AC=1;AF=0.50;AN=2;BaseQRankSum=3.241;DP=40;Dels=0.00;FS=0.000;HRun=1;HaplotypeScore=0.0000;MQ=50.61;MQ0=1;MQRankSum=-4.823;QD=2.65;ReadPosRankSum=2.054;SB=-66.05;snp128    GT:AD:DP:GQ:PL  0/1:27,13:40:99:136,0,718
# sz2317@b1a.c1.titan 120509_SN828_0132_AC0J91ACXX$ head -2 mapping/1.bam_refine/VarCalling/list.vcf-files.txt.vcf.dbsnp.annovar.exonic_variant_function
# line1   synonymous SNV  Rrs1:NM_021511:exon1:c.T630G:p.T210T,   chr1    9536234 9536234 T       G       chr1    9536234 rs31477593      T       G       26.78   LowQual AB=0.730;AC=1;AF=0.50;AN=2;BaseQRankSum=-1.453;DP=37;Dels=0.00;FS=9.688;HRun=2;HaplotypeScore=0.9967;MQ=52.19;MQ0=0;MQRankSum=-4.497;QD=0.72;ReadPosRankSum=-0.804;SB=-0.01;snp128   GT:AD:DP:GQ:PL  0/1:27,10:37:56.77:57,0,768
# line2   synonymous SNV  Rrs1:NM_021511:exon1:c.A648G:p.E216E,   chr1    9536252 9536252 A       G       chr1    9536252 rs13467800      A       G       105.92  PASS    AB=0.675;AC=1;AF=0.50;AN=2;BaseQRankSum=3.241;DP=40;Dels=0.00;FS=0.000;HRun=1;HaplotypeScore=0.0000;MQ=50.61;MQ0=1;MQRankSum=-4.823;QD=2.65;ReadPosRankSum=2.054;SB=-66.05;snp128    GT:AD:DP:GQ:PL  0/1:27,13:40:99:136,0,718
#
#
my %hash = ();
for ($count = 1; $count <= 19; $count++) {
 	$hash{ "chr$count" } = ();
 }
       $hash{ "chrX" } = ();
       $hash{ "chrY" } = ();
       $hash{ "chrM" } = ();

open(fin1, "<" . $fin1) || die("Could not open $fin1 file!");
my @field=();
while(<fin1>) {
	chomp;
	@field = split ("\t", $_);
	$hash{ $field[2] }{ $field[3] } = [ $field[7], $field[8], $field[9], $field[10], $field[11], $field[12], $field[13], $field[14].";refGene.function=$field[0];refGene.geneName=$field[1]", $field[15], $field[16]];
}
close fin1;

# print  "\n".join("\t", @{$hash{"chr1"}{"9536234"}} );

open(fin2, "<" . $fin2) || die("Could not open $fin2 file!");
my @field=();
while(<fin2>) {
        chomp;
        @field = split ("\t", $_);
	if (exists $hash{ $field[3] }{ $field[4] }){
#	print "\n $field[3]	$field[4]	$field[1]	$field[2] ";
		$field[2]  =~ s/,$//g;
        	$hash{ $field[3] }{ $field[4] }[7].=";refGene.functionalClass=$field[1];refGene.name=$field[2]";
	}	
	else {
		print "Warning:Does not exist $field[3]  $field[4] !!! \n";
		exit;
	}
}
close fin2;

open(fin3, "<" . $fin3) || die("Could not open $fin3 file!");
while(<fin3>) {
	if ($_ =~ m/^#/) {  print $_; }
	else { last; }
}
close fin3;


my $firstime=0;
for ($count = 1; $count <= 19; $count++) {
	foreach $key (sort {$a <=> $b} keys %{$hash {"chr$count"}}) {
		if ($firstime > 0) {print "\n";}
		print join("\t", @{$hash{"chr$count"}{$key}});
		$firstime = $firstime + 1;
	}
}

foreach $key (sort {$a <=> $b} keys %{$hash {"chrM"}}) {
	print "\n".join("\t", @{$hash{"chrM"}{$key}});
}

foreach $key (sort {$a <=> $b} keys %{$hash {"chrX"}}) {
		print "\n".join("\t", @{$hash{"chrX"}{$key}});             
}

foreach $key (sort {$a <=> $b} keys %{$hash {"chrY"}}) {
                print "\n".join("\t", @{$hash{"chrY"}{$key}});
}

