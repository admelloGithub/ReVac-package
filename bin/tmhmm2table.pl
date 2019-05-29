#!/usr/bin/perl -w
use strict;

my $tmhmm_out = $ARGV[0];    #tmhmm raw output list of files
my $tmhmm2table	= $ARGV[1];   #Output file for 2 column table
my @files;

if ($tmhmm_out =~ /.list/) {
	open (my $list, "<", $tmhmm_out) || die "Couldn't open tmhmm raw list file:$!\n";
	@files = <$list>;
	close $list || die "$!\n";
} else {
	die "TMHMM inputs should be a list of files\n";
}

if (-e $tmhmm2table) {
    system ("rm $tmhmm2table");
    system ("touch $tmhmm2table");
}

foreach (@files) {
	
	my @pred_tmhs;
	my $tmhmm_raw = $_;
	$tmhmm_raw =~ s/\s+$//;
	#print "$tmhmm_raw\n";

	open (my $fh, "<", $tmhmm_raw) || die "Couldn't open tmhmm raw output file:$!\n";

	while (<$fh>) {
		if ($_ =~ /^# / && $_ =~ /Number of predicted TMHs:/) {
			push (@pred_tmhs, $_);
		}
	}

	close $fh || die "$!\n";

	open (my $fh2, ">>", $tmhmm2table) || die "Couldn't open tmhmm2table output file:$!\n";

	#print $fh2 "Sequence_ID\tNo_of_TMHs\n";

	foreach (@pred_tmhs) {
    	my @split = split (/\s+/, $_);
    	my @seqid_tag = split (/\./, $split[1]);
		print $fh2 "$seqid_tag[0]\t$split[6]\n";
	}	

	close $fh2 || die "$!\n";

}
