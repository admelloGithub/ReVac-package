#!/usr/bin/perl -w
use strict;
use List::MoreUtils qw(uniq);

my $netctlpan_cutoff = $ARGV[0];
my @immuno_input;

open (my $fh, "<", $netctlpan_cutoff) || die "$!\n";
while (<$fh>) {
	my @line = split(/\t/, $_);
	push (@immuno_input, $line[1]);
}
close $fh || die "$!\n";

my @input=uniq(@immuno_input);

foreach (@input) {
print "$_\n";
}
