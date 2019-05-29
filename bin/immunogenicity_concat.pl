#!/usr/bin/perl -w
use strict;
use lib "/usr/local/projects/ergatis/package-revac/lib/modules";
use List::MoreUtils qw(uniq);

if (scalar(@ARGV) > 2 || scalar(@ARGV) == 0) {
	die "Requires 2 files, an MHC_I immuno input and a NetCTLpan immuno input, the order is irrelevant.\n";
}

my @peptides;

foreach (@ARGV) {
	open (my $fh, "<", $_) || die "Couldnt open :$_ because \"$!\"\n";
		@peptides = <$fh>;
	close $fh || die "$!\n";
}

my @uniq = uniq(@peptides);

foreach (@uniq) {
	print "$_\n" if ($_ !~ /[X]+/) ;
}
