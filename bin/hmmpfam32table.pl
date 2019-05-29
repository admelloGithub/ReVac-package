#!/usr/bin/perl -w
use strict;

my $hmmpfam_out = $ARGV[0]; #Must be a list file of raw outputs
my $hmmpfam32table = $ARGV[1]; #Output file for 2 table column
my (@files,%table);

if ($hmmpfam_out =~ /.list$/) {
    open (my $list, "<", $hmmpfam_out) || die "Couldn't open HMMPFAM3 raw list file :$!\n";
    @files = <$list>;
    close $list || die "$!\n";
} else {
    die "HMMPFAM3 output should be a list file \".list\" of raw outputs\n";
}

if (-e $hmmpfam32table) {
	system ("rm $hmmpfam32table");
    system ("touch $hmmpfam32table");
}

foreach (@files) {

	my $hmmpfam_raw = $_;
	$hmmpfam_raw =~ s/\s+$//;
	my $text;

	open (my $fh, "<", $hmmpfam_raw) || die "Couldn't open hmmpfam raw output file: $!\n";

	read $fh, $text, -s $fh;

	close $fh || die "$!\n";

	my @pep_hmms = split (/\/\//, $text);

	foreach (@pep_hmms) {

    	my @lines = split (/\n/, $_);

			for (my $i=0; $i<scalar(@lines); $i++) {

	    		if ($lines[$i] =~ /^Query:/) {
  				my @idline = split (/\s+/, $lines[$i]);
				my @evalues = split (/[ ]{2,}/, $lines[$i+5]);
					if ($evalues[1] || $evalues[10]) {
		    			$evalues[10] =~ s/:$//;
		    			$table{$idline[1]} = "$evalues[1]|$evalues[10]";
					} else {
		    			$table{$idline[1]} = "None";
					}
	    		}

			}	

	}

}

open (my $fh2, ">>", $hmmpfam32table) || die "Couldn't open hmmpfam32table output file: $!\n";

#print $fh2 "Sequence_ID\tE-value/Protein\n";

foreach (keys %table) {
	my @split = split(/\./, $_);
	print $fh2 "$split[0]\t$table{$_}\n";
}

close $fh2 || die "$!\n";


