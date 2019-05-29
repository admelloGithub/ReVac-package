#!/usr/bin/perl -w
use strict;

use File::Basename;

my $list = $ARGV[0];
my $outdir = $ARGV[1];
my %contigs;
open (my $fh, "<", $list) || die "$!\n";
        my @files = <$fh>;
close $fh || die "$!\n";

foreach (@files) {
	$_ =~ s/\s+$//;
	open (my $fh2, "<", $_) || die "File:$_\t$!\n";
		my $id;
		while (<$fh2>) {
			if ($_ =~ />/) { #gnl\|.*\|(.*)\s(.*)/) {
			$id = $_;
			$id =~ s/^.*\|/>/;
			#$id =~ s/\s+/_/g;
			#print $id."\n";
			} else {
			$contigs{$id} .= "$_";
			}
		}
	close $fh2;
		
}

foreach (sort keys %contigs) {
	my $file = $_;
	$file =~ s/ .*$//;
	$file =~ s/>//;
	open (my $fh3, ">", "$outdir/$file.fasta") ||die "$!\n";
		print $fh3 "$_\n$contigs{$_}";
	close $fh3;
}
