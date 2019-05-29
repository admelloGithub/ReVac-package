#!/usr/bin/perl -w
use strict;

my $seq_input = $ARGV[1];  #Must be a multifasta file of all input sequences used in phobos
my $phobos_out = $ARGV[0]; #Must be a list file of phobos outputs
my $phobos2table = $ARGV[2]; #Output file for 2 column table
my (@files,%table);

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

if ($phobos_out =~ /.list/) {
    open (my $list, "<", $phobos_out) || die "Couldn't open phobos raw list :$!\n";
    @files =<$list>;
    close $list || die "$!\n";
} else {
    die "Phobos file must be a list file \".list\" of raw outputs \n";
}

foreach (@files) {

	my $phobos_raw = $_;
	$phobos_raw =~ s/\s+$//;
	my @lines;

	open (my $fh2, "<", $phobos_raw) || die "Couldn't open phobos raw output file:$!\n";

	while (<$fh2>) {

    	if ($_ !~ /^#/) {
			push (@lines, $_);
    	}

	}

	close $fh2 || die "$!\n";

	shift @lines;

	foreach (@lines) {

    	my @split = split(/\s+/, $_);
    	$split[0] =~ s/^\s+//;
    	$split[0] =~ s/\s+$//;
    	$table{$split[0]} = "$split[2]/$split[6]/$split[12]";

	}

}

if (-e $phobos2table) {
	system ("rm $phobos2table");
    system ("touch $phobos2table");
}

open (my $fh3, ">", $phobos2table) || die "Couldn't open phobos2table output file:$!\n";

#print $fh3 "Sequence_ID\tPerfection/Repeat_No/Repeat_unit\n";

foreach (keys %table) {
	my @split = split(/\./, $_);
	print $fh3 "$split[0]\t$table{$_}\n";
}

close $fh3 || die "$!\n";

