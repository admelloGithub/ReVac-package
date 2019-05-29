#!/usr/bin/perl -s

# $Id: residues,v 1.6 1997/03/31 15:28:08 lixinz Exp lixinz $

if ($h) {
    print "residues - number of nt/aa of each fasta sequence\n";
    print "Usage: residues multi_seq_fasta_db\n";
    exit (0);
}

$/="\n>";
while (<>) {
    chop if />$/;
    s/>// if /^>/;  ### $_ = ">" . "$_" if ! /^>/;
    ($gene,$seq) = /\s*(\S*).*\n([\s\S]*)/;
    $seq =~ s/\s+//g;
##    $tot = $seq =~ tr/a-zA-Z//;
    $tot = $seq =~ tr/a-zA-Z\*//;
    $grandtot += $tot;
    $t = $seq =~ tr/tT//;
    $c = $seq =~ tr/cC//;
    $g = $seq =~ tr/gG//;
    $a = $seq =~ tr/aA//;
    $n = $seq =~ tr/aAcCgGtT//c;
    $p = (($c+$g)/$tot)*100;
    if ($n) {
	printf "$gene\t$tot\t$grandtot\tperG+C:%2.1f\tT:$t\tC:$c\tG:$g\tA:$a\tother:$n\n",$p;
    } else {
	printf "$gene\t$tot\t$grandtot\tperG+C:%2.1f\tT:$t\tC:$c\tG:$g\tA:$a\n",$p;
    }
}

__END__
