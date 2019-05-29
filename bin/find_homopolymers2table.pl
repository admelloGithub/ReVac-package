#!/usr/bin/perl -w
use strict;

my $find_homops_out = $ARGV[0];  #Must be a list file of raw outputs
my $find_homops2table = $ARGV[1];	#Output file for 2 table column
my (@files,%table);

if ($find_homops_out =~ /.list$/) {
	open (my $list, "<", $find_homops_out) || die "Couldn't open find_homopolymers raw list file :$!\n";
	@files = <$list>;
	close $list || die "$!\n";
} else {
	die "Find_homopolymers output should be a list file \".list\" of raw outputs\n";
}

foreach (@files) {

	my $find_homops_raw = $_;
	$find_homops_raw =~ s/\s+$//;

	open (my $fh, "<", $find_homops_raw) || die "Couldn't open find_homopolymers raw output file :$!\n";

	my $text;

	while (<$fh>) {
    	$text = do { local $/; <$fh> };
	}

	close $fh || die "$!\n";

	my @seqs = split(/>/, $text);

	foreach (@seqs) {
    	my @homops;
    	my @split = split(/\s+/, $_);
		for (my $i=4; $i<scalar(@split); $i+=3) { 
	    	my @nucleotide = split(//, $split[$i]);
	    	my $l = scalar(@nucleotide); 
	    	push(@homops, "$nucleotide[0]$l");
	    	$table{$split[0]} = join("|",@homops);
		}
	}

}

if (-e $find_homops2table) {
    system ("rm $find_homops2table");
    system ("touch $find_homops2table");
}

open (my $fh2, ">>", $find_homops2table) || die "Couldn't open find_homopolymers2table output file :$!\n";
#	print $fh2 "Sequence_ID\tNo_of_homopolymers\n";

foreach (keys %table) {
	my @split = split(/\./, $_);
	print $fh2 "$split[0]\t$table{$_}\n";
}

close $fh2 || die "$!\n";

