#!/usr/bin/perl -w
use strict;

my @newlines;
my $file = $ARGV[0];
my %cogs;

open (my $fh, "<", $file) || die "$!\n";

my @lines = <$fh>;

close $fh || die "$!\n";

for (my $i = 1; $i<scalar(@lines); $i++) {

	my @split = split(/\t/,$lines[$i]);
	my $id = shift @split;
	$id =~ s/\.CDS\.[0-9]+\.1//;
	foreach (@split) {
		if ($_ >= 0.60) {
			$_ = 1;
			$cogs{$id} += 1;
		} else {
			$_ = 0;
		}
	}

	my $nline = $id."\t".join("\t",@split);
	push (@newlines,$nline);	

}
		
#print "$lines[0]";
#foreach (@newlines) {
#	print "$_\n";
#}

foreach (keys %cogs) {
	print "$_\t$cogs{$_}\n";
}

