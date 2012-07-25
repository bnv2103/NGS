#!/usr/bin/perl

#For Mouse : Convert gene annotation from Annovar back to VCF 
my $fin1 = $ARGV[0];   # variant_function file
my $fin2 = $ARGV[1];   # exonic.variant_function file
my $fin3 = $ARGV[2];   # original VF file, to grab header 

#VCF files assuming chr 1 to 22 and X and Y only
#headerssss
#chr1       3421959 other fields..
#chr1       3688030 other fields..

# head -2 mapping/1.bam_refine/VarCalling/list.vcf-files.txt.vcf.dbsnp.annovar.variant_function
# exonic  Rrs1    chr1    9536234 9536234 T       G       chr1    9536234 rs31477593      T       G       26.78   LowQual AB=0.730;AC=1;AF=0.50;AN=2;BaseQRankSum=-1.453;DP=37;Dels=0.00;FS=9.688;HRun=2;HaplotypeScore=0.9967;MQ=52.19;MQ0=0;MQRankSum=-4.497;QD=0.72;ReadPosRankSum=-0.804;SB=-0.01;snp128   GT:AD:DP:GQ:PL  0/1:27,10:37:56.77:57,0,768
# exonic  Rrs1    chr1    9536252 9536252 A       G       chr1    9536252 rs13467800      A       G       105.92  PASS    AB=0.675;AC=1;AF=0.50;AN=2;BaseQRankSum=3.241;DP=40;Dels=0.00;FS=0.000;HRun=1;HaplotypeScore=0.0000;MQ=50.61;MQ0=1;MQRankSum=-4.823;QD=2.65;ReadPosRankSum=2.054;SB=-66.05;snp128    GT:AD:DP:GQ:PL  0/1:27,13:40:99:136,0,718
# sz2317@b1a.c1.titan 120509_SN828_0132_AC0J91ACXX$ head -2 mapping/1.bam_refine/VarCalling/list.vcf-files.txt.vcf.dbsnp.annovar.exonic_variant_function
# line1   synonymous SNV  Rrs1:NM_021511:exon1:c.T630G:p.T210T,   chr1    9536234 9536234 T       G       chr1    9536234 rs31477593      T       G       26.78   LowQual AB=0.730;AC=1;AF=0.50;AN=2;BaseQRankSum=-1.453;DP=37;Dels=0.00;FS=9.688;HRun=2;HaplotypeScore=0.9967;MQ=52.19;MQ0=0;MQRankSum=-4.497;QD=0.72;ReadPosRankSum=-0.804;SB=-0.01;snp128   GT:AD:DP:GQ:PL  0/1:27,10:37:56.77:57,0,768
# line2   synonymous SNV  Rrs1:NM_021511:exon1:c.A648G:p.E216E,   chr1    9536252 9536252 A       G       chr1    9536252 rs13467800      A       G       105.92  PASS    AB=0.675;AC=1;AF=0.50;AN=2;BaseQRankSum=3.241;DP=40;Dels=0.00;FS=0.000;HRun=1;HaplotypeScore=0.0000;MQ=50.61;MQ0=1;MQRankSum=-4.823;QD=2.65;ReadPosRankSum=2.054;SB=-66.05;snp128    GT:AD:DP:GQ:PL  0/1:27,13:40:99:136,0,718
#
#
my %hash = ();
my @otherVariants = ();  ##holds other variants on same locus as in %hash
# push @otherVariants, ["null","value","here"];

for ($count = 1; $count <= 19; $count++) {
 	$hash{ "chr$count" } = ();
 }
       $hash{ "chrX" } = ();
       $hash{ "chrY" } = ();
       $hash{ "chrM" } = ();

## Proces Annovar.variant_function file
open(fin1, "<" . $fin1) || die("Could not open $fin1 file!");
my @field=();
while(<fin1>) {
	chomp;
	@field = split ("\t", $_);
	## chromosome and pos from Annovar.variant_function file (at col 7 ,8)
	#first check if this is not already present in $hash then add as a new entry 
#	print "$field[7]\t$field[8]\t";
	if (not exists $hash{ $field[7] }{ $field[8] } )
	{
		$hash{ $field[7] }{ $field[8] } = [ $field[7], $field[8], $field[9], $field[10], $field[11], $field[12], $field[13], $field[14].";refGene.function=$field[0];refGene.geneName=$field[1]", $field[15], $field[16]];
#		print "ZERO\n";
	}
	elsif  (join(" ",@{$hash{ $field[7] }{ $field[8] }}[0..6]) ne join(" ",@{[ $field[7], $field[8], $field[9], $field[10], $field[11], $field[12], $field[13]]} ))
	{ ## now cheick if is not simply a duplicate entry , only then add a new entry in @othervVariants
#		print STDERR "$field[7]\t$field[8]\t",join("_",@{$hash{ $field[7] }{ $field[8] }}[0..6]), "\t", join("_",@{[ $field[7], $field[8], $field[9], $field[10], $field[11], $field[12], $field[13]]} ), "\n";
#		print "ZERO";	
		my @tempArr=@{$hash{ $field[7] }{ $field[8] }};
		my $flag=0;		
		foreach $index (@{$tempArr[10..(scalar(@tempArr)-1)]}) {
#			print "chk$index";
			if (join(" ",@{$otherVariants[$index]}[0..6]) eq join(" ",@{[ $field[7], $field[8], $field[9], $field[10], $field[11], $field[12], $field[13]]} ))
			{
				$flag=1; 
#				print "fail$index"; 
				last;
			}
		}
		if ($flag == 0){## only make new entry if a duplicate entry was not found in @otherVariants
			push @otherVariants, [ $field[7], $field[8], $field[9], $field[10], $field[11], $field[12], $field[13], $field[14].";refGene.function=$field[0];refGene.geneName=$field[1]", $field[15], $field[16]];
			push @{$hash{$field[7]}{$field[8]}}, scalar(@otherVariants)-1 ;
#			print "pass", scalar(@otherVariants) , "\n";
		}else{
#                       print "Warning:Duplicate\n"; # entry for $field[7]  $field[8] of $fin1 !!! \n";
		}
	}else{
#	$hash{ $field[7] }{ $field[8] } = [ $field[7], $field[8], $field[9], $field[10], $field[11], $field[12], $field[13], $field[14].";refGene.function=$field[0];refGene.geneName=$field[1]", $field[15], $field[16], 0];
#		print "ZERO:Warning:Duplicate\n"; #entry for $field[7]  $field[8] of $fin1 !!! \n";
#		exit;	
	}
}

close fin1;

#foreach $item (@otherVariants)
#{
#	print STDERR join("\t", @{$item})."\n";
#}
# exit;
## Process Annovar.exonic_variant_function file
open(fin2, "<" . $fin2) || die("Could not open $fin2 file!");
my @field=();
while(<fin2>) {
        chomp;
        @field = split ("\t", $_);
	if (exists $hash{ $field[8] }{ $field[9] }){	## chromosome and pos from Annovar.exonic_variant_function file (at col 8,9 )
		## search for the entry in %hash
		if  (join(" ",@{$hash{ $field[8] }{ $field[9] }}[0..6]) eq join(" ",@{[  $field[8], $field[9], $field[10], $field[11], $field[12], $field[13], $field[14]]} ))
		{
			$field[2]  =~ s/,$//g;
        		$hash{ $field[8] }{ $field[9] }[7].=";refGene.functionalClass=$field[1];refGene.name=$field[2]";
		}
		else ##if cannot find in %hash, then search if it is in anotherVariants
		{
			my @tempArr=@{$hash{ $field[8] }{ $field[9] }};
                	foreach $index (@{$tempArr[10..(scalar(@tempArr)-1)]}) {
                        if (join(" ",@{$otherVariants[$index]}[0..6]) eq join(" ",@{[ $field[8], $field[9], $field[10], $field[11], $field[12], $field[13], $field[14]]} ))
                        {	##if you find the variant, update its refGene.functionalClass etc and break .
				$otherVariants[$index][7].=";refGene.functionalClass=$field[1];refGene.name=$field[2]";
                                last;
                        }}
                }
	}	
	else {
		print "Warning:Does not exist $field[8]  $field[9] at $field[0] of $fin2 !!! \n";
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
		print join("\t", @{$hash{"chr$count"}{$key}}[0..9]);
		my $itemref = \@{$hash{"chr$count"}{$key}};
		foreach $index (@{$itemref}[10..(scalar(@{$itemref})-1)]) {
#			print STDERR $index."\n";
			print "\n".join("\t",@{$otherVariants[$index]});
		}
		$firstime = $firstime + 1;
	}
}

foreach $key (sort {$a <=> $b} keys %{$hash {"chrM"}}) {
	print "\n".join("\t", @{$hash{"chrM"}{$key}}[0..9]);
		my $itemref = \@{$hash{"chrM"}{$key}};
		foreach $index (@{$itemref}[10..(scalar(@{$itemref})-1)]) {
			print "\n".join("\t",@{$otherVariants[$index]});
		}
}

foreach $key (sort {$a <=> $b} keys %{$hash {"chrX"}}) {
		print "\n".join("\t", @{$hash{"chrX"}{$key}}[0..9]);             
		my $itemref = \@{$hash{"chrX"}{$key}};
		foreach $index (@{$itemref}[10..(scalar(@{$itemref})-1)]) {
			print "\n".join("\t",@{$otherVariants[$index]});
		}
}

foreach $key (sort {$a <=> $b} keys %{$hash {"chrY"}}) {
                print "\n".join("\t", @{$hash{"chrY"}{$key}}[0..9]);
		my $itemref = \@{$hash{"chrY"}{$key}};
		foreach $index (@{$itemref}[10..(scalar(@{$itemref})-1)]) {
			print "\n".join("\t",@{$otherVariants[$index]});
		}
}

