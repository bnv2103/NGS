#!/usr/bin/env perl
#
use warnings;
use strict;

while (<>) {
	chomp;
	my @parts = split /\t/;
	if ( $parts[10] == 1 ){
	print "@","$parts[0]:$parts[2]:$parts[3]:$parts[4]:$parts[5]#$parts[6]/$parts[7]\n";
	print "$parts[8]\n";
	print "+","$parts[0]:$parts[2]:$parts[3]:$parts[4]:$parts[5]#$parts[6]/$parts[7]\n";
	print "$parts[9]\n";
	}
	else{
        warn "@","$parts[0]:$parts[2]:$parts[3]:$parts[4]:$parts[5]#$parts[6]/$parts[7]\n";
        warn "$parts[8]\n";
        warn "+","$parts[0]:$parts[2]:$parts[3]:$parts[4]:$parts[5]#$parts[6]/$parts[7]\n";
        warn "$parts[9]\n";
        }
}
