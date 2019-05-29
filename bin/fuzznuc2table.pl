#!/usr/bin/perl -w
use strict;

my $seq_input = $ARGV[1];   #Multifasta file of all sequences used in fuzznuc 
my $fuzznuc_out = $ARGV[0];	#List file of all fuzznuc raw outputs
my $fuzznuc2table = $ARGV[2];	#Output file for 2 column table
my (%table,@files);

open (my $fh, "<", $seq_input) || die "Couldn't open input sequence file :$!\n";

while (<$fh>) {
    if ($_ =~ /^>/) {
	my $strip = $_;
	$strip =~ s/>//g;
	$strip =~ s/^\s+//;
	$strip =~ s/\s+$//;
	$table{$strip} = "None";
    }
}

close $fh || die "$!\n";

if ($fuzznuc_out =~ /.list/) {
	open (my $list, "<", $fuzznuc_out) || die "Couldn't open fuzznuc raw list :$!\n";
	@files =<$list>;
	close $list || die "$!\n";
} else {
	die "Fuzznuc file must be a list file \".list\" of raw outputs \n";
}

foreach (@files) {

	my $fuzznuc_raw = $_;
	$fuzznuc_raw =~ s/\s+$//;

	open (my $fh2, "<", $fuzznuc_raw) || die "Couldn't open fuzznuc raw output file :$!\n";

	while (<$fh2>) {
    	if ($_ =~ /^Name: /) {
        	my $strip = $_;
        	$strip =~ s/^Name: //g;
        	$strip =~ s/^\s+//;
        	$strip =~ s/\s+$//;
        	$table{$strip} = "Match";
    	}
	}

	close $fh2 || die "$!\n";

}

if (-e $fuzznuc2table) {
    system ("rm $fuzznuc2table");
    system ("touch $fuzznuc2table");
}

open (my $fh3, ">>", $fuzznuc2table) || die "Couldn't open fuzznuc2table output file :$!\n";

#print $fh3 "Sequence_ID\tFuzznuc\n";

foreach (keys %table) {
	my @split = split(/\./, $_);
	print $fh3 "$split[0]\t$table{$_}\n";
}

close $fh3 || die "$!\n";

