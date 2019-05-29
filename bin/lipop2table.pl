#!/usr/bin/perl -w
use strict;

my $lipop_out = $ARGV[0];  #List of all raw output
my $lipop2table = $ARGV[1]; #Output filepath for 2 column table
my (@files,%table);

if ($lipop_out =~ /.list/) {
    open (my $list, "<", $lipop_out) || die "Couldn't open LipoP raw output list for reading:$!\n";
    @files = <$list>;
    close $list || die "$!\n";
} else {
    die "LipoP output should be a list file \".list\"\n";
}

if (-e $lipop2table) {
    system ("rm $lipop2table");
    system ("touch $lipop2table");
}

foreach (@files) {

	my $lipop_raw = $_;
	$lipop_raw =~ s/\s+$//;

	open (my $fh, "<", $lipop_raw) || die "Couldn't open Lipop raw output file:$!\n";

	while (<$fh>) {
    	my @line = split(/\s+/, $_);
    	$line[3] =~ s/score=//g;
    	$table{$line[1]} = "$line[2]|$line[3]";
	}

	close $fh || die "$!\n";

}

open (my $fh2, ">>", $lipop2table) || die "Couldn't open Lipop2table file for writing:$!\n";

foreach (keys %table) {
	my @split = split(/\./, $_);
    print $fh2 "$split[0]\t$table{$_}\n";
}

close $fh2 || die "$!\n";
