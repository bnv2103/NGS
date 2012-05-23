#!/usr/bin/perl

#Sprt 1 VCF file
my $fin1 = $ARGV[0]; #infile with chr & pos in 1st two columns
my $sift = $ARGV[1];  #sift score file
my $pp2  = $ARGV[2];  #pp2 file

open(fin1, "<" . $fin1) || die("Could not open $fin1  file!");

my %hash = ();
for ($count = 1; $count <= 22; $count++) {
 	$hash{ $count } = ();
 }
       $hash{ 'X' } = ();
       $hash{ 'Y' } = ();

my @field=();
#Parse chromo file fin1
while(<fin1>) {
  if ($_ =~ m/^#/) {
  	next;
  }
  else{
	chomp;
	@field = split ('\t', $_);
	$field[0] =~ s/chr//i; 
       #  print "$field[0]\t $field[1]\n";

	if ( $filed[0]  =~ m/\d/  or  $filed[0]  =~ m/x/i or $filed[0] =~ m/y/i ) {}else{
		$hash{ $field[0] }{ $field[1] } = ['',''];
#		print "$field[0]\t $field[1]\n";
	}
  }
}
close fin1;

# print  %hash ;
## PArse SIFT file
#Chr     Pos     REF     ALT     SIFT_score
#1       100174670       C       T       0.33
open(fin1, "<" .  $sift) || die("Could not open $sift SIFT  file!");

while(<fin1>) {
        chomp;
        @field = split ('\t', $_);
 	if ( $field[0] =~ m/chr/i ) {next; }
	if ( defined $hash{ $field[0] }{ $field[1] } ) {
		$hash{ $field[0] }{ $field[1] }[0] =  $field[4] ;
        }
}
close fin1;

## PArse PP2 file
#Chr     Pos     REF     ALT     PP2_score
#1       865628  G       A       0.955
open(fin1, "<" .  $pp2) || die("Could not open $pp2 PP2 file!");

while(<fin1>) {
        chomp;
        @field = split ('\t', $_);
 	if ( $field[0] =~ m/chr/i ) {next; }
	if ( defined $hash{ $field[0] }{ $field[1] } ) {
		$hash{ $field[0] }{ $field[1] }[1] = $field[4] ;
        }
}
close fin1;


#Print in sorted order 
print "Chromosome\tPosition\tSIFT\tPP2\n";
for ($count = 1; $count <= 22; $count++) {
	foreach $key (sort {$a <=> $b} keys %{$hash {$count}}) {
		print "$count\t$key\t$hash{$count}{$key}[0]\t$hash{$count}{$key}[1]\n";
	}
}

foreach $key (sort {$a <=> $b} keys %{$hash {'X'}}) {
                print "X\t$key\t$hash{'X'}{$key}[0]\t$hash{'X'}{$key}[1]\n";
}

foreach $key (sort {$a <=> $b} keys %{$hash {'Y'}}) {
                print "Y\t$key\t$hash{'Y'}{$key}[0]\t$hash{'Y'}{$key}[1]\n";
}

