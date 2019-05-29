#!/usr/bin/perl -w
use strict;

my $blat_raw_list = $ARGV[0];
my $outdir = $ARGV[1];
open (my $FH, "<", $blat_raw_list) || die "$!\n";
my @files = <$FH>;
close $FH || die "$!\n";

foreach (@files) {
	my (%seqs,%matches,%regions,$id);
	my $blat_file = $_;
	$blat_file =~ s/\s+//;
	
	open (my $fh,"<",$blat_file) || die "$!:$blat_file\n";

	while (<$fh>) {
	  
		if ($_ =~ /^Query=/) {
			$id = $_;
			$id =~ s/^Query= //;
			$id =~ s/\n//;
			$id =~ s/\.polypeptide\.[0-9]+\.1//;
		} elsif ($_ =~ /^Query:/) {
			$seqs{$id} .= "$_|";
		}

	}

	close $fh || die "$!\n";

	foreach my $id (keys %seqs) {
		my @split = split(/\|/,$seqs{$id});
			foreach (@split) {
				my @split = split(/\s+/,$_);
				if ($matches{$id}{"$split[1]"}) {
					my @oldsplit = split(/\t/,$matches{$id}{"$split[1]"});
					if ($oldsplit[1] < $split[3]) {
					$matches{$id}{"$split[1]"} = "$split[1]\t$split[3]\t$split[2]";
					}
				} else {
					$matches{$id}{"$split[1]"} = "$split[1]\t$split[3]\t$split[2]";
				}
			}
	}

	foreach my $id (keys %matches) {

		my @regions;
		my $peptide = "";	

		foreach (sort {$a<=>$b} keys %{$matches{$id}}) {
			#print "$id\t$matches{$id}{$_}\n";
			push (@regions, $matches{$id}{$_});
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
                        my $str = substr($split2[2],$split1[1]-$split2[0]+1);
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
	open (my $fh2, ">>", "$outdir/blat2table.out") || die "$!\n";
	foreach (keys %regions) {
		my @split = split(/\./,$regions{$_});
		pop @split;
		my $l = length(join("",@split));
		print $fh2 "$_\t".scalar(@split)."($l)\n";
	}
	close $fh2 || die "$!\n";
}
