#!/usr/bin/perl -w
use strict;

my $netctlpan_exec = $ARGV[0];
my $method = $ARGV[1];
my $length = $ARGV[2];
my $threshold = $ARGV[3];
my $cleavage = $ARGV[4];
my $tap = $ARGV[5];
my $epitope = $ARGV[6];
my $i_file_path = $ARGV[7];
my @cmds;

my @alleles = qw(HLA-A01:01 HLA-A02:01 HLA-A03:01 HLA-A24:02 HLA-A26:01 HLA-B07:02 HLA-B08:01 HLA-B27:05 HLA-B39:01 HLA-B40:01 HLA-B58:01 HLA-B15:01);

foreach (@alleles) {
	push (@cmds, "$netctlpan_exec\t$method\t$length\t$_\t$threshold\t$cleavage\t$tap\t$epitope\t$i_file_path");
}

foreach (@cmds){
	print "$_\n";
}






