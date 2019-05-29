#!/usr/bin/perl -w
use strict;

my $psortb_out = $ARGV[0]; #Should be a list of psortb outputs in long format 
my $psortb2table = $ARGV[1]; #Output file for 2 column table
my (@files,%table);

if ($psortb_out =~ /.list/) {
	open (my $list, "<", $psortb_out) || die "Couldn't open psortb raw output list for reading:$!\n";
   	@files = <$list>;
    close $list || die "$!\n";
} else {
    die "PSORTB output should be a list file \".list\"\n";
}

foreach (@files) {

	my $psortb_raw = $_;
	$psortb_raw =~ s/\s+$//;

	open (my $fh, "<", $psortb_raw) || die "Couldn't open psortb raw output file :$!\n";

	my $headers = <$fh>;
	$headers =~ s/_Score//g;
	my @headers = split(/\t/, $headers);

	#if (scalar(@headers) != 35 || scalar(@headers) != 30) {
    #	die "Potentially incorrect psortb format. Psortb raw output should be in long format.\n";
	#}

	if (scalar(@headers) == 35) {		###Gram negative output

		while (<$fh>) {
			my @line = split(/\t/, $_);
			my %scores;
			if ($line[30] =~ /^Unknown/) {
				for (my $i=25; $i<30; $i++) {                             ##This is to report highest scores. Not final prediction.
					if ($scores{$line[$i]}) {
						my $var = ($i/10000);								##Handles duplicate keys (same score) removed later.
						$scores{$line[$i]+$var} = $headers[$i];
					} else {
						$scores{$line[$i]} = $headers[$i];
					}
					if ($line[$i] >= 7.5) {
						$table{$line[0]}="$line[$i]|$headers[$i]";
					}
				}
				my @scores = sort keys %scores;
				if (!($table{$line[0]})) {
					if ($scores[$#scores-1] + $scores[$#scores] >= 7.5) {
                        my $s= sprintf("%.2f",$scores[$#scores]);
                        my $sc = sprintf("%.2f",$scores[$#scores-1]);
                        $table{$line[0]} ="$sc|$scores{$scores[$#scores-1]},$s|$scores{$scores[$#scores]}";
                    } elsif ($scores[$#scores-1] + $scores[$#scores] + $scores[$#scores-2]> 7.6) {
                        my $s= sprintf("%.2f",$scores[$#scores]);
                        my $sc = sprintf("%.2f",$scores[$#scores-1]);
                        my $sco = sprintf("%.2f",$scores[$#scores-2]);
                        $table{$line[0]} ="$sco|$scores{$scores[$#scores-2]},$sc|$scores{$scores[$#scores-1]},$s|$scores{$scores[$#scores]}";
                    }
				}
			}
			if ($line[30] !~ /^Unknown/) {
				$table{$line[0]} = $line[30];
			}

		}

	} elsif (scalar(@headers) ==30) {

		while (<$fh>) {
            my @line = split(/\t/, $_);
            my %scores;
            if ($line[25] =~ /^Unknown/) {
                for (my $i=21; $i<25; $i++) {                             ##This is to report highest scores. Not final prediction.
                    if ($scores{$line[$i]}) {
						my $var = ($i/10000);
                        $scores{$line[$i]+$var} = $headers[$i];
                    } else {
	                    $scores{$line[$i]} = $headers[$i];
                    }
                    if ($line[$i] >= 7.5) {
                        $table{$line[0]}="$line[$i]|$headers[$i]";
                    }
                }
                my @scores = sort keys %scores;
                if (!($table{$line[0]})) {
                    if ($scores[$#scores-1] + $scores[$#scores] >= 7.5) {
						my $s= sprintf("%.2f",$scores[$#scores]);
						my $sc = sprintf("%.2f",$scores[$#scores-1]);
                        $table{$line[0]} ="$sc|$scores{$scores[$#scores-1]},$s|$scores{$scores[$#scores]}";
                    } elsif ($scores[$#scores-1] + $scores[$#scores] + $scores[$#scores-2]> 7.6) {
						my $s= sprintf("%.2f",$scores[$#scores]);
                        my $sc = sprintf("%.2f",$scores[$#scores-1]);
                        my $sco = sprintf("%.2f",$scores[$#scores-2]);
                        $table{$line[0]} ="$sco|$scores{$scores[$#scores-2]},$sc|$scores{$scores[$#scores-1]},$s|$scores{$scores[$#scores]}";
                    }
                }
             }
             if ($line[25] !~ /^Unknown/) {
                 $table{$line[0]} = $line[25];
             }

         }

	} else {
		die "Potentially incorrect psortb format. Psortb raw output should be in long format.\n";
	}

	close $fh || die "$!\n";

}

if (-e $psortb2table) {
    system ("rm $psortb2table");
    system ("touch $psortb2table");
}

open (my $fh2, ">", $psortb2table) || die "Couldn't open psortb2table file for writing:$!\n";

foreach (keys %table) {
	my @split = split(/\./, $_);
    print $fh2 "$split[0]\t$table{$_}\n";
}

close $fh2 || die "$!\n";
