#!/usr/bin/perl -w
use strict;

#use File::Slurp;

#Takes antigenic raw output and generates 2 tab delimited columns of seq id \t "YES. Score=XXXX / NO"

my $antigenic_output = $ARGV[0];
my $antigenic2table = $ARGV[1];
my @seqid_lines;
my ($file, @files);

if ($antigenic_output =~ /.list$/) {
	open (my $list, "<", $antigenic_output) || die "Couldn't open antigenic raw output list: $!\n";
	@files = <$list>;
	close $list || die "$!\n";
} else {
	die "Antigenic output must be a list of raw files\n";
} 

if (-e $antigenic2table) {
    system ("rm $antigenic2table");
	system ("touch $antigenic2table");
}

foreach (@files) {

	my $antigenic_raw = $_;
	$antigenic_raw =~ s/\s+$//;
	#print "$antigenic_raw\n";
	open (my $fh, "<", $antigenic_raw) || die "Couldn't open antigenic raw output file: $!\n";

	read $fh, $file, -s $fh;
	#$file = read_file($fh);

	close $fh;

	my @split = split(/#---------------------------------------\n#---------------------------------------/, $file);

	foreach (@split) {
	
   		my @lines = split(/\n/, $_);
	
		foreach (@lines) {

	   		if ($_ =~ /^# Sequence:/) {
				push (@seqid_lines, $_);
	 	    }

			elsif ($_ =~ /^Sequence:/) {
                my $score = $_;
                $score =~ s/Sequence: //g;
                push (@seqid_lines, $score);
            }

	   		elsif ($_ =~ /^Score:/) {
				my $score = $_;
				$score =~ s/Score: //g;
				push (@seqid_lines, $score);
	   		}

		}

	}


	my $string = join (' ',@seqid_lines);

	my @id_n_scores = split(/#/, $string);
	shift @id_n_scores;
	
	if (-e $antigenic2table) {
		system ("rm $antigenic2table");
		system ("touch $antigenic2table");
	}

	open (my $fh2, ">>", $antigenic2table) || die "Couldn't open antigenic2table output file: $!\n";

	#print $fh2 "Sequence_ID\tAntigenicity\n";

	foreach (@id_n_scores) {

    	my @final = split (/\s+/, $_);

		if ($final[2] && $final[7]) {	
    	    my @split = split(/\./, $final[2]);
			my $l = scalar(@final);
			$l = ($l-7)/2;
			my $coverage;
				for(my $i=7; $i<scalar(@final); $i+=2) {			
					$coverage .= $final[$i];
				}
			my $cov_len = length($coverage);
			print $fh2 "$split[0]\t$final[7]|$final[8]|$l($cov_len)\n";
		} 
		elsif ($final[2]) {
			my @split = split(/\./, $final[2]);
	    	print $fh2 "$split[0]\tNone\n";
		}
		else { print "Antigenic Format may be different from when this script was written.\n"; }

	}

	close $fh2;

}
