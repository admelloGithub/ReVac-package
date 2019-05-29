#!/usr/bin/perl -w
use strict;
use File::Basename;

my $bcell_all_list = $ARGV[0]; 
my $bcell2table = $ARGV[1];
my (%seqs,$seq);

open (my $FH, "<", $bcell_all_list) || die "$!\n";
        my @bcellfiles = <$FH>;
close $FH || die "$!\n";

open (my $fh2, ">", $bcell2table) || die "$!\n";

foreach (@bcellfiles) {

	my $bcell_all = $_;
	$bcell_all =~ s/\s+$//;
	my $lines;
	open (my $fh, "<", $bcell_all) || die "$bcell_all.$!\n";

	while (<$fh>) {
		$lines .= $_;
	}
	my @split = split(/>/,$lines);
	shift @split;			

	close $fh || die "$!\n";

	foreach (@split) {
		my @split1 = split(/\n/,$_);
		$seqs{$split1[0]} = join("\n",@split1[1..$#split1]);
	}

}

foreach my $id (keys %seqs) {
	my @regions = split(/\n/,$seqs{$id});
	shift @regions;	
	my $peptides = "";
	my $peptide = "";

		for(my $i=0; $i<=scalar(@regions)-1; $i++) {

			if ($i == $#regions) {
			   my @split =split(/\s+/,$regions[$i]);
			   if ($peptide) { $peptides .= "$peptide."; } else {
					$peptides .= "$split[$#split].";
			   }
			} else {
				my $line1 = $regions[$i];
				my $line2 = $regions[$i+1];
				my @split1 = split(/\s+/,$line1);
				my @split2 = split(/\s+/,$line2);
				$peptide = $split1[$#split1] if ($peptide eq "");
					if (!$peptide) {print STDERR "$line1\n";}
						if ($split1[1]>=$split2[0]) {
							 my $str = substr($split2[$#split2],$split1[1]-$split2[0]+1);
							 $peptide .= $str;
						} else {
							 $peptides .= "$peptide.";
							 $peptide = "";
						}

			}

		}

	my @peptides = split(/\./,$peptides);
	my $length;
	my $ct=1;
	my $locus_tag = $id;
    $locus_tag =~ s/^>//;
    $locus_tag =~ s/\.polypeptide\.[0-9]+\.1//;
		foreach (@peptides) {
			$length += length($_);
			print ">$locus_tag\_$ct\n$_\n";		##Prints out regions in a multifasta file
			$ct++;
		}
	print $fh2 "$locus_tag\t".scalar(@peptides)."($length)|".scalar(@regions)."\n" if (scalar(@regions) != 0);
	#print $fh2 join (".",@peptides)."\n" if (scalar(@regions) != 0);

}

close $fh2 || die "$!\n";
