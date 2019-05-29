#!/usr/local/bin/perl -w
use strict;
use List::MoreUtils qw(uniq);

my $mhc_i_cutoff = $ARGV[0];
my $immuno_out = $ARGV[1];
my $immuno2table = $ARGV[2];
my $immuno_peptides= $ARGV[1];
$immuno_peptides =~ s/\.raw$/\.peptides/;
my (%seqs,%table,%peptides,%alleles,%count,%regions,%imm_peps,%hash);

open (my $fh, "<", $mhc_i_cutoff) || die "Couldn't open mhc class I cutoff file for reading:$!\n";
#    my @lines = <$fh>;
#close $fh || die "Couldnt close:$!\n";

#foreach (@lines) {
while (<$fh>) {
	my @split = split(/\s+/,$_);
	$split[0] =~ s/\.p.*$//;			####To bypass netctlpan id trimming
	$split[1] =~ s/^\s+//; $split[1] =~ s/\s+$//;
	$seqs{$split[1]} .= "$split[0]."; 
	$regions{$split[0]}{$split[1]}{start} = $split[3];
	$regions{$split[0]}{$split[1]}{end} = $split[4];
	$alleles{$split[1]} .= "$split[5].";
}
close $fh || die "$!\n";

foreach (keys %alleles) {
	my @split = split(/\./, $alleles{$_});
	my @uniq = uniq @split;
	#pop @uniq;
	$alleles{$_} = join(".",@uniq);
}

foreach (keys %seqs) {
    my @split = split(/\./, $seqs{$_});
    my @uniq = uniq @split;
	#pop @uniq;
	$seqs{$_} = join(".",@uniq);
}

#foreach my $id (keys %regions) {
#	foreach my $pep (keys %seqs) {
		#push (@imm_peps, "$id\t$pep\t$regions{$id}{$pep}{start}\t$regions{$id}{$pep}{end}");
#	foreach my $id (keys %regions) {
#		$imm_peps{$id}{$regions{$id}{$pep}{start}} = "$id\t$pep\t$regions{$id}{$pep}{start}\t$regions{$id}{$pep}{end}";
#	}
#	}
#}
#foreach (keys %seqs) {
#print $_.";$seqs{$_};\n";
#}


open (my $fh2, "<", $immuno_out) || die "Couldn't open immuno out file for reading $immuno_out :$!\n";

while (<$fh2>) {
    if ($. > 4) {
	my @cols = split(/,/,$_);
		if ($cols[2] >= 0.2) {	
	   		if ($seqs{$cols[0]}) {
				my @seqids = split(/\./,$seqs{$cols[0]});
					foreach (@seqids) {
						if ($_ eq "") {next;}
						$table{$_} += 1;
						$count{$_} .= "$alleles{$cols[0]}.";
						$peptides{$_} .= "$cols[0].";
						$imm_peps{$_}{$regions{$_}{$cols[0]}{start}} = "$_\t$cols[0]\t$regions{$_}{$cols[0]}{start}\t$regions{$_}{$cols[0]}{end}";
					}
#		   		if ($table{$seqs{$cols[0]}}) {
#				$table{$seqs{$cols[0]}} += 1; 	
#				} else {
#				$table{$seqs{$cols[0]}} += 1;
#				}
#				$count{$seqs{$cols[0]}} .= $alleles{$cols[0]}; 
	    	} else { die "Output peptide isn't in cutoff file :$cols[0]:"; }
#	    	$peptides{$seqs{$cols[0]}} .= "$cols[0].";		
    	}
    }
}
close $fh2 || die "$!\n";

foreach my $id (keys %imm_peps) {
    #foreach my $allele (keys %{$hash{$id}}) {
      my @regions;
      my $peptide = "";

        foreach my $start (sort {$a <=> $b} keys %{$imm_peps{$id}} ) {
            push (@regions, $imm_peps{$id}{$start});
        }

        for(my $i=0; $i<=scalar(@regions)-1; $i++) {            ###Compares sorted peptide start and end postions to determine if overlapping for each protein and
                                                                ###concatenates them unitl the next one doesn't overlap.
            if ($i == $#regions) {
                my @split =split(/\s+/,$regions[$i]);
                if ($peptide) { $hash{$split[0]} .= "$peptide.|"; } else {
                    $hash{$split[0]} .= "$split[1].|";
                }
            } else {

                my $line1 = $regions[$i];
                my $line2 = $regions[$i+1];
                my @split1 = split(/\s+/,$line1);
                my @split2 = split(/\s+/,$line2);
                $peptide = $split1[1] if ($peptide eq "");
                if ($split1[3]>=$split2[2]) {
                    my $str = substr($split2[1],$split1[3]-$split2[2]+1);
                    $peptide .= $str;
                } else {
                    $hash{$split1[0]} .= "$peptide.";
                    $peptide = "";
                }

            }

        }

    #}

}


open (my $fh3, ">", $immuno2table) || die "Couldn't open immunogenicity2table file for writing:$!\n";

foreach (keys %table) {
	my @split = split(/\./, $count{$_});
	my @uniq = uniq @split;
	my $lengtha = scalar(@uniq);
	my $alleles = join(".",@uniq);
	my @splitr = split(/\./,$hash{$_});
	pop @splitr;
	my $reg_cov = length ( join ("",@splitr) );
	my $lengthr = scalar(@splitr);

	print $fh3 "$_\t$lengthr($reg_cov)|$table{$_}|$lengtha\n";			#|$alleles\n";
}

close $fh3 ||die "$!\n";

open (my $fh4, ">", $immuno_peptides) || die "Couldn't open immunogenic peptides file for writing:$!\n";

foreach (keys %peptides) {
	print $fh4 "$_\t$peptides{$_}\nRegions_$_:\t$hash{$_}\n";
}

close $fh4 ||die "$!\n";
