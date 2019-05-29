#!/usr/local/bin/perl -w
use strict;

my $map_file=$ARGV[0];
my $cog_file=$ARGV[1];
my %hash;

open (my $fh, "<", $map_file) || die "$!\n";

my @ids = <$fh>;

close $fh;

my %ids;
foreach (@ids) {
	my @split = split(/\t/,$_);
	$split[1] =~ s/^\s+//;
    $split[1] =~ s/\s+$//;
	$ids{$split[1]}=$split[0];
}

open (my $fh2, "<", $cog_file) || die "$!\n";
my $cogs;
while (<$fh2>) {
$cogs .= $_;
}
close $fh2;

my @cogs = split(/COG = /,$cogs);
shift @cogs;

foreach (@cogs) {
	my @lines = split(/\n/,$_);
	$lines[0] =~ s/,.*//;
	my $id = "j_ortholog_cluster_$lines[0]";
	my $ct = scalar(@lines) - 1;
	shift @lines;
	foreach (@lines) {
		$_ =~ s/^\s+//;
		$_ =~ s/\s+$//;
		#print "$_\n";
		$hash{$ids{$_}} = "$id|$ct";
	}
}

foreach (keys %hash) {
	print "$_\t$hash{$_}\n";
}
