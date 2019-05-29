#!/usr/bin/perl -w
use strict;

my $list = $ARGV[0];
my %lengths;

open (my $FILE,"<", $list) or die "Could not open $list:$!\n";
my @files = <$FILE>;
close $FILE;

foreach my $f (@files) {

	$f =~ s/\s+$//;
	open (my $fh,"<", $f) or die "Could not open $f:$!\n";
		my @lines = <$fh>;
	close $fh;
	
	foreach (@lines) {
		my @split = split(/\s+/,$_);
		$lengths{$split[0]} = $split[1];
	}

}

foreach (sort keys %lengths) {
	print "$_\t$lengths{$_}\n";
}
