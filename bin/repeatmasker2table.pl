#!/usr/bin/perl -w
use strict;

my $repeat_masked_fsa = $ARGV[0];  #Must be a list of masked sequences from Repeatmasker
#my $repeat_raw = $ARGV[1];
my $repeat2table = $ARGV[1]; #Output for 2 column table
my (%table,@files);

if ($repeat_masked_fsa =~ /.list/) {
	open (my $list, "<", $repeat_masked_fsa) || die "Couldn't open repeatmasker raw output list for reading:$!\n";
    @files = <$list>;
    close $list || die "$!\n";
} else {
    die "Repeatmasker output should be a list file \".list\"\n";
}

if (-e $repeat2table) {
    system ("rm $repeat2table");
    system ("touch $repeat2table");
}

foreach (@files) {

	my ($text,%seqs);
	my $repeat_masked_raw = $_;
	$repeat_masked_raw =~ s/\s+$//;
	
	if (-e $repeat_masked_raw) {
    	open (my $fh, "<", $repeat_masked_raw) || die "Couldn't open repeatmasker masked file:$!\n";
		$text = do { local $/; <$fh> };
    	close $fh || die "$!\n";
	} else {
    	exit;
	}

	my @seqs = split (/>/, $text);
	shift @seqs;

	foreach (@seqs) {
    	my @split = split (/\n/, $_);
    	my $seq = join('', @split[1..(scalar(@split)-1)]);
		if ($seq =~ m/N+/) {                               ##Confirms only masked sequences are matched to their repeats
    	    $seqs{$split[0]} = $seq; 
		}
	}

	my $repeat_raw = $repeat_masked_raw;
	$repeat_raw =~ s/masked$/out/;

	open (my $fh2, "<", $repeat_raw) || die "Couldn't open repeatmasker raw outfile for reading:$!\n";

	my $line = <$fh2>;
	if ($line =~ /no repetitive sequences/) {
    	next;
	}

	while (<$fh2>) {
    	if ($. > 3) {
		my @split = split(/\s+/, $_); 
	    	if (defined $seqs{$split[5]}) {
	    		$table{$split[5]} .= "|$split[10]|";
	    	}
    	}
	}

	close $fh2 || die "$!\n";

}

open (my $fh3, ">>", $repeat2table) || die "Couldn't open repeatmasker2table file for writing:$!\n";

foreach (keys %table) {
	my @split = split(/\./, $_);
	print $fh3 "$split[0]\t$table{$_}\n";
}

close $fh3 || die "$!\n";
