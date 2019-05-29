#!/usr/bin/perl -w
use strict;

my $homologs_input = $ARGV[0];
my $autoimmunity2table = $ARGV[1];
my (%seqs,%coverage,%regions,$text);

open (my $fh, "<", $homologs_input) || die "Couldn't open autoimmunity homologous peptides file :$!\n";

while (<$fh>) {
	$text .= $_;
}

close $fh || die "$!\n";

my @seqs = split (/>/, $text);
shift @seqs;

foreach (@seqs) {
	my @lines = split(/\n/,$_);
	my @ids = split(/,/,$lines[2]);
		foreach my $line (@lines) {
			if ($line =~ /\^[GALMFWKQESPVICYHRNDTX]+\^/) {
				#$line =~ s/^[0-9]+\s+\^//;
				#$line =~ s/\^$//;
				foreach my $id (@ids) {
					$seqs{$id} += 1;
					$coverage{$id} .= "$line.";
				}
			}
		}
}

foreach my $id (keys %coverage) {

	my @peps = split(/\./,$coverage{$id});
	my %coord;
	foreach (@peps) {
		my @split = split(/\s+/,$_);
		$split[1] =~ s/[0-9]+\s+\^//;
		$split[1] =~ s/\^$//;
		$split[1] =~ s/^\^//;
		$split[0]++;
		my $l = length($split[1])+$split[0];
		my $key = "$split[0].$l";					###Same peptides are overwritten in the hash. Key is as decimal start.end to allow sorting later.
		$coord{$key} = "$split[0]\t$l\t$split[1]";
	}
	
	 my @regions;
	 my $peptide = "";

        foreach my $start (sort {$a <=> $b} keys %coord ) {
            push (@regions, $coord{$start});
        }

        for(my $i=0; $i<=scalar(@regions)-1; $i++) {            ###Compares sorted peptide start and end postions to determine if overlapping for each protein and
                                                                ###concatenates them unitl the next one doesn't overlap.    
            if ($i == $#regions) {
                my @split =split(/\s+/,$regions[$i]);
                if ($peptide) { $regions{$id} .= "$peptide.|"; } else {
                    $regions{$id} .= "$split[2].|";
                }
            } else {
 
                my $line1 = $regions[$i];
                my $line2 = $regions[$i+1];
                my @split1 = split(/\s+/,$line1);
                my @split2 = split(/\s+/,$line2);
                $peptide = $split1[2] if ($peptide eq "");
                if ($split1[1]>=$split2[0]) {
                    if (($split1[1]-$split2[0])<=($split2[1]-$split2[0])) {
						my $str = substr($split2[2],$split1[1]-$split2[0]);       
                    	$peptide .= $str;
					} else {
						$peptide =~ s/$split2[2].+$/$split2[2]/;
						next;
					}
                } else {
                    $regions{$id} .= "$peptide.";
                    $peptide = "";
                }

            }

        }

}

foreach my $key (keys %seqs) {
	my @split = split(/\./,$key);
	$split[0] =~ s/^\s+//;
	my $regions = $regions{$key};
	my @regions = split(/\./,$regions);
	pop @regions;
	my $no_of_reg = scalar(@regions);
	$regions =~ s/[\|\.]//g;
	my $l = length($regions);
	print "$split[0]\t$no_of_reg($l)\n";#$regions{$key}\n";
}
