#!/usr/bin/perl -w
use strict;
use List::MoreUtils qw(uniq);

my $immuno_input = $ARGV[0]; ###Takes file of a column of peptides and ensures there's only 1 copy of each.

open (my $fh, "<", $immuno_input) || die "$!\n";

my @peps = <$fh>;

close $fh || die "$!\n";

my @uniq = uniq @peps;

open (my $fh2, ">", $immuno_input) || die "$!\n";

print $fh2 @uniq;
	
close $fh2 || die "$!\n";
