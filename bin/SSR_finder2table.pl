#! /usr/bin/perl -w
use strict;

my $ssr_genes_file = $ARGV[0];
my %ssrs;

open (my $fh, "<", $ssr_genes_file) || die "$!\n";

while (<$fh>) {
	my @split = split(/\s+/,$_);
	my @motif = split(/_/,$split[0]);
	$ssrs{"N$split[4]"} .= "/$split[3]|$split[1]|$motif[0]($motif[1])/";
}
close $fh || die "$!\n";

foreach (sort keys %ssrs) {
print "$_\t$ssrs{$_}\n";
}
