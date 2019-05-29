#! /usr/bin/perl -w

use strict;
use File::Basename;

my $blastfile = $ARGV[0];
my $outdir = $ARGV[1];
my $pepfile = $ARGV[2];
my @idsnseqs;

my $basename = basename($blastfile);


open (my $fh, "<", $blastfile) || die "Couldn't open blast file:$!\n";

my @lines = <$fh>;

close $fh;

foreach (@lines) {

	if ($_ =~ /^Query=/ || $_ =~ /^>/) {
		push (@idsnseqs, $_);
	}

}

open (my $fh3, "<", $pepfile) || die "Couldn't open Homologous Peptides File :$!\n";

my $peplines;
while(<$fh3>){
        $peplines.=$_;
}
my @peplines = split(/^>/m, $peplines);
shift @peplines;

close $fh3;

unless (-e "$outdir/$basename.out") {
`touch $outdir/$basename.out`;
}

open (my $fh2, ">", "$outdir/$basename.out") || die "Couldn't open out file :$!\n";

for (my $i=0; $i <= scalar(@idsnseqs); $i++) {
    if ($idsnseqs[$i] =~ m/^Query=/ and $idsnseqs[$i+1] =~ m/^>/) {
	print $fh2 "$idsnseqs[$i]";
    } elsif ($idsnseqs[$i] =~ /^>/ and $idsnseqs[$i+1] =~ m/^Query=/) {
        print $fh2 "\t$idsnseqs[$i]\n";
    } elsif ($idsnseqs[$i] =~ /^>/) {
	print $fh2 "\t$idsnseqs[$i]";
    }
}

foreach (@peplines) {
print $fh2 ">$_\n";
}
close $fh2;













