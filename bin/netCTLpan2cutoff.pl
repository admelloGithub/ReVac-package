#!/usr/bin/perl -w
use strict;
use File::Basename;
my $netctlpan_list = $ARGV[0];
my %peptides;

open (my $fh, "<", $netctlpan_list) || die "$!\n";
my @files = <$fh>;
close $fh || die "$!\n";

my $outdir = dirname($netctlpan_list);

foreach my $file (@files) {

	my $allele = basename($file);
	$allele =~ s/\.netctlpan_all\.out\s+//;
	$allele =~ s/.*\.//;
	$file =~ s/\s+$//;
	open (my $fh2, "<", $file) || die "$!\n";

	my $id;
	while (<$fh2>) {

		if ($_ =~ /:$/) {
			$id = $_;
			$id =~ s/:\s+//;
		} elsif ($_ =~ /^[0-9]+/) {
			my @split = split(/\s+/,$_);
			$peptides{$id}{$split[0]}{$allele} = "$split[1]" if ($split[$#split] <=1.0);
		}	

	}

	close $fh2 || die "$!\n";

}

open (my $fh3, ">", "$outdir/Netctlpan") || die "$!\n";
foreach my $id (keys %peptides) {
    foreach my $co (sort keys %{$peptides{$id}}) {
        foreach my $alle (keys %{$peptides{$id}{$co}}) {
			my $end = $co+(length($peptides{$id}{$co}{$alle})-1);
            print $fh3 "$id\t$peptides{$id}{$co}{$alle}\t-\t$co\t$end\t$alle\n";
        }
    }
}
close $fh3 || die "$!\n";
