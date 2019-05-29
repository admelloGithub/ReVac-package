#!/usr/bin/perl -w
use strict;

my $input = $ARGV[0];
my $output = $ARGV[1];
my @peptides;
my @ranks;
my %hash;

open (my $fh, "<", $input) || die "Couldn't open input mhc_class_ii_all.out file:$!\n";

while (<$fh>) {
	if ($. > 2) {
	my @split = split(/\s+/, $_);
	$hash{$split[5]} = $split[6];
	}
}

close $fh;

open (my $fh2, ">", $output) || die "Couldn't open output file:$!\n";

foreach (keys %hash) {
    if ($hash{$_} <= "5.00") {
	print $fh2 "$_\n";
    }
}

close $fh2;

