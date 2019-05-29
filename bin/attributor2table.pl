#!/usr/bin/perl -w
use strict;

if (scalar(@ARGV) == 0) {
	print "No command line attributor output files entered. \nPlease provide any number of files as cmd line arguments.\n";
}

foreach(@ARGV){

my %seqs;
my $file = $_;

open (my $fh, "<", $file) || die "$!\n";

while (<$fh>) {
	if ($_ =~ /^>/) {
		my $fhead = $_;
		my @split = split(/\s+/,$fhead);
		$split[0] =~ s/^>//;
		if ( $split[$#split] =~ /^go::/ ){
			my $GO = pop @split;
			$GO =~ s/^go:://;
			$seqs{$split[0]} = join("_",@split[1..$#split])."|$GO";
		} else {
			$seqs{$split[0]} = join("_",@split[1..$#split]);
		}
	}
}

close $fh || die "$!\n";

foreach (sort keys %seqs) {
	print "$_\t$seqs{$_}\n";
}

}
