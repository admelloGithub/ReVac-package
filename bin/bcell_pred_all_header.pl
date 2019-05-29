#!/usr/bin/perl -w
use strict;

my $file = $ARGV[0];

open (IN, $file) || die "Cannot open file ".$file." for read";     
my @lines=<IN>;  
close IN;
 
open (OUT, ">", $file) || die "Cannot open file ".$file." for write";
for (my $i=0; $i<scalar(@lines); $i++)
{  
   #my $j= $i+1;
   print OUT "$lines[$i]";  
}
close OUT;


system ("sed -i '1i \$\;I_FILE_BASE\$\;\t\$\;PYTHON\$\;\t\$\;BCELL\$\;\t\$\;METHOD\$\;\t\$\;FILE_PATH\$\;' $file");
