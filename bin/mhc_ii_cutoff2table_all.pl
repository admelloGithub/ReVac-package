#!/usr/bin/perl -w
use strict;
use List::MoreUtils qw(uniq);

my $mhc_cutoff = $ARGV[0];
my (%hash,%allele,%regions,%peptides);
my $mhc_cutoff2table = $mhc_cutoff."2table.out";

open (my $fh, "<", $mhc_cutoff) || die "$!\n";

while (<$fh>) {
	my @split = split(/\s+/, $_);
	my $allele = $split[5];
		if ($allele{$split[0]}) { 	
			$allele{$split[0]} .= "$split[5]."; 		#Stores all alleles per protein, including repeated ones.
			$peptides{$split[0]} .= "$split[1].";		#Stores all peptides per protein, including repeated ones.
		} else {
			$allele{$split[0]} = "$split[5].";
			$peptides{$split[0]} = "$split[1].";
		}
	$hash{$split[0]}{$split[3]} = $_;  					#Stored as hash{SeqID}{Peptide_start_pos}=Line 
}														#For each SeqID the Peptide_start_pos will overlap for other alleles, to store only unique peptides in the hash

close $fh || die "$!\n";

foreach my $id (keys %hash) {
	#foreach my $allele (keys %{$hash{$id}}) {
	  my @regions;
	  my $peptide = "";

		foreach my $start (sort {$a <=> $b} keys %{$hash{$id}} ) {
			push (@regions, $hash{$id}{$start});
		}
	
		for(my $i=0; $i<=scalar(@regions)-1; $i++) {			###Compares sorted peptide start and end postions to determine if overlapping for each protein and
																###concatenates them unitl the next one doesn't overlap.	
			if ($i == $#regions) {
				my @split =split(/\s+/,$regions[$i]);
				if ($peptide) { $regions{$split[0]} .= "$peptide.|"; } else {
					$regions{$split[0]} .= "$split[1].|";
				}
			} else {
	
				my $line1 = $regions[$i];
				my $line2 = $regions[$i+1];
				my @split1 = split(/\s+/,$line1);
				my @split2 = split(/\s+/,$line2);
				$peptide = $split1[1] if ($peptide eq "");
				if ($split1[4]>=$split2[3]) {
					my $str = substr($split2[1],$split1[4]-$split2[3]+1);		
					$peptide .= $str;
				} else {
					$regions{$split1[0]} .= "$peptide.";
					$peptide = "";
				}

			}

		}

	#}	
	
}
													###Actual sequence data is present in the hashes(regions,peptides & alleles). Print them to access actual data.

open (my $fh2, ">", $mhc_cutoff2table) || die "$!\n";		

foreach (keys %regions) {
	my @splita = split(/\./,$allele{$_});
	my @splitp = split(/\./,$peptides{$_});
	my @splitr = split(/\./,$regions{$_});
	pop @splitr;
	my @uniqa = uniq @splita;
	my @uniqp = uniq @splitp;
	my $lengtha = scalar(@uniqa);
	my $lengthp = scalar(@uniqp);
	my $lengthr = scalar(@splitr);
    my $l;
		foreach (@splitr) {
        	$l += length($_);
        }
	print $fh2 "$_\t$lengthr($l)|$lengthp|$lengtha\n";
}

close $fh2 || die "$!\n";

#EOF
