#!/usr/bin/perl

## Function: position phv repeats with regards to ORFs: in promoter, 5' third, middle, or 3' third.

if ($#ARGV != 1) {
		print "\nNeeds 2 parameters: phv.list orf.list\n\n"
		    . "phv.list: feat_name end5 end3\n"
		    . "orf.list: feat_name end5 end3 com_name\n\n";
		exit;
		}

($phvs,$orfs) = @ARGV;

open (PHVS, "$phvs") || die "Cannot open $phvs : $!\n";
while(<PHVS>) {
	chomp;
	($phv,$pend5,$pend3,$seq) = split;
	$coord = (($pend3-$pend5)/2)+$pend5;
	open (ORFS, "$orfs") || die "Cannot open $orfs : $!\n";
	while(<ORFS>) {
	    chomp;
	    ($orf,$oend5,$oend3,$com_name) = split(/\t/);
	    if ($oend3>$oend5) {
		if ($coord>=$oend5-50 && $coord<=$oend3) {
		    $third=($oend3-$oend5)/3;
		    if ($coord>=$oend5-50 && $coord<$oend5) {
			print "$phv\t$pend5\t$pend3\tpromot\t$orf\t$oend5\t$oend3\t$com_name\t$seq\n";
		    }
		    elsif ($coord>=$oend5 && $coord<$oend5+$third) {
			print "$phv\t$pend5\t$pend3\t5prime\t$orf\t$oend5\t$oend3\t$com_name\t$seq\n";
		    }
		    elsif ($coord>=$oend5+$third && $coord<$oend5+(2*$third)) {
			print "$phv\t$pend5\t$pend3\tmiddle\t$orf\t$oend5\t$oend3\t$com_name\t$seq\n";
		    }
		    else {
			print "$phv\t$pend5\t$pend3\t3prime\t$orf\t$oend5\t$oend3\t$com_name\t$seq\n";
		    }
		}
	    }
	    else {
		if ($coord>=$oend3 && $coord<=$oend5+50) {
		    $third=($oend5-$oend3)/3;
		    if ($coord<=$oend5+50 && $coord>$oend5) {
			print "$phv\t$pend5\t$pend3\tpromot\t$orf\t$oend5\t$oend3\t$com_name\t$seq\n";
		    }
		    elsif ($coord<=$oend5 && $coord>(2*$third)+$oend3) {
			print "$phv\t$pend5\t$pend3\t5prime\t$orf\t$oend5\t$oend3\t$com_name\t$seq\n";
		    }
		    elsif ($coord<=(2*$third)+$oend3 && $coord>$third+$oend3) {
			print "$phv\t$pend5\t$pend3\tmiddle\t$orf\t$oend5\t$oend3\t$com_name\t$seq\n";
		    }
		    else {
			print "$phv\t$pend5\t$pend3\t3prime\t$orf\t$oend5\t$oend3\t$com_name\t$seq\n";
		    }
		}
	    }
	}
	close (ORFS) || die "Cannot close $orfs : $!\n";
    }
close (PHVS) || die "Cannot close $phvs : $!\n";

__END__

