#!/usr/bin/perl -w
use strict;

my $epitope_out = $ARGV[0];
my $immuno_in = $ARGV[1];
my %epi_pep;

open (my $fh, "<", $epitope_out) || die "Couldn't open epitope_cluster raw output :$!\n";

while (<$fh>) {

    my @split = split(/\s+/,$_);
    $epi_pep{$split[3]} = $split[4];

}

close $fh;

foreach (keys %epi_pep) {
	print "$_\t$epi_pep{$_}\n";
}
