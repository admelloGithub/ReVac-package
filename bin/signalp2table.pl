#!/usr/bin/perl -w
use strict;

#Takes signalp raw output and generates 2 tab delimited columns of seq id \t Signal_peptide amino acid sequence
my $signalp_output = $ARGV[0];   #List of all raw output
my $input_fasta = $ARGV[1];		#Single multifasta file of all sequences used to generate signalp outputs
my $signalp2table = $ARGV[2];	#Output filepath for 2 column table
my (@files,%seqs,%positions,%table);

if ($signalp_output =~ /.list/) {
	open (my $list, "<", $signalp_output) || die "Couldn't open signalp raw output list for reading:$!\n";
	@files = <$list>;
	close $list || die "$!\n";
} else {
	die "Signalp output should be a list file \".list\"\n";
}

if (-e $signalp2table) {
	system ("rm $signalp2table");
	system ("touch $signalp2table");
}

open (my $fh, "<", $input_fasta) || die "Couldn't open input multifasta file for reading:$!\n";
    my $text = do { local $/; <$fh> };
close $fh || die "$!\n";

my @seqs = split (/>/, $text);
shift @seqs;

foreach (@seqs) {
    my @split = split (/\n/, $_);
    $seqs{$split[0]} = join('', @split[1..(scalar(@split)-1)]);
}

foreach (@files) {

	my $signalp_raw = $_;
	$signalp_raw =~ s/\s+$//;
	
	if ($signalp_raw =~ /.gz$/) { 
		system ("gunzip $signalp_raw");
		$signalp_raw =~ s/.gz$//;
	}

	#print "$signalp_raw\n";

	open (my $fh2, "<", $signalp_raw) || die "Couldn't open signalp raw output file: $!\n";

    my $lines = do { local $/; <$fh2> };
    my @split = split(/#/,$lines); 
	foreach (@split) {
	   if ($_ =~ /YES/) {
		my @sp = split(/\n/,$_);
		my @array = split(/\s+/,$sp[5]);
		my $position = $array[2] if ($array[5] =~ /YES/);
		my @array2 = split(/\s+/,$sp[6]);
		my $seq_id = $array2[0];
		$seq_id =~ s/^Name=//;
		$positions{$seq_id} = $position if ($array[5] =~ /YES/);
	   }
	}

	close $fh2 || die "Couldn't close signalp raw output file :$!\n";

	foreach (keys %positions) {
    	if ($seqs{$_}) {
		my @position = split(/-/,$positions{$_});
		my $sp = substr($seqs{$_}, $position[0]-1, $position[1]);
		$table{$_} = $sp;
    	}
	}
	system ("gzip $signalp_raw");

}

open (my $fh3, ">>", $signalp2table) || die "Couldn't open signalp2table.out file:$!\n";

foreach (keys %table) {
	my @split = split(/\./, $_);
	print  $fh3 "$split[0]\t$table{$_}\n";
}

close $fh3 || die "$!\n";
