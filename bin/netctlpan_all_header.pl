#!/usr/bin/perl -w
use strict;

my $file = $ARGV[0];

open (IN, $file) || die "Cannot open file ".$file." for read";     
my @lines=<IN>;  
close IN;
 
open (OUT, ">", $file) || die "Cannot open file ".$file." for write";
for (my $i=0; $i<scalar(@lines); $i++)
{  
   my $j= $i+1;
   print OUT "$j\t$lines[$i]";  
}
close OUT;


system ("sed -i '1i \$\;I_FILE_BASE\$\;\t\$\;CTLPAN_EXEC\$\;\t\$\;METHOD\$\;\t\$\;LENGTH\$\;\t\$\;ALLELE\$\;\t\$\;THRESHOLD\$\;\t\$\;TAP_WEIGHT\$\;\t\$\;CLEAVAGE_WEIGHT\$\;\t\$\;EPITOPE_THRESHOLD\$\;\t\$\;FILE_PATH\$\;' $file");
