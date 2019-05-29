#!/usr/bin/perl -w
use strict;

my $matchtable = $ARGV[0];
my %clusters;

open (my $fh, "<", $matchtable) || die "$!\n";
while (<$fh>) {
	my @line = split(/\t/,$_);
	my $c_size;
	my $cid = $line[0];
	shift @line;
	foreach my $lt (@line) {
        if ($lt =~ /---/) {next;}
		else {
			$c_size += 1;
		}
	}
	foreach my $id (@line) {
		if ($id =~ /---/) {next;}
		elsif ($clusters{$id}) {
			$id =~ s/\s+$//;
			$clusters{$id} .= ",PanOCT_cluster_$cid|$c_size";
		} else {
			$id =~ s/\s+$//;
			$clusters{$id} = "PanOCT_cluster_$cid|$c_size";
		}
	}
}
close $fh || die "$!\n";

foreach (sort keys %clusters) {
	print "$_\t$clusters{$_}\n";
} 
