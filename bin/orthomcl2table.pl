#!/usr/bin/perl -w
use strict;

my $orthomcl_groups = $ARGV[0];
my %cogs;
open (my $fh, "<", $orthomcl_groups) || die "$!\n";

while (<$fh>) {
	my @split = split(/\s+/,$_);
	if (scalar(@split) <= 2) { next; }
	$split[0] =~ s/://;
	my $cluster_id = $split[0];
	shift @split;
	foreach my $lt (@split) {
		$cogs{$lt} = $cluster_id."|".scalar(@split);
	}
}

close $fh;

foreach (sort keys %cogs) {
	print "$_\t$cogs{$_}\n";
}
