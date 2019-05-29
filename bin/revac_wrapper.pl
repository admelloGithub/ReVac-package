#!/usr/bin/perl -w
use strict;
use Getopt::Std;

getopts("o:a:n:c:t:d:");
our ($opt_o,$opt_a,$opt_n,$opt_c,$opt_t,$opt_d);
# -o yes/no
# -a yes/no
# -n number_of_genomes
# -c conservation_ratio(0.8)
# -t summary_table_file
# -d out_dir
`/usr/local/projects/ergatis/package-revac/bin/revac_algorithm.pl $opt_t $opt_a $opt_n > $opt_d/Summary_table.scored.txt`;

if ($opt_o eq "yes") {

`/usr/local/projects/ergatis/package-revac/bin/revac_recall_orthologs_withssr.pl $opt_d/Summary_table.scored.txt > $opt_d/Summary_table.scored.recalled.txt`;

system("for i in `cat $opt_d/Summary_table.scored.recalled.txt|grep -v \"#\"|cut -f3|sort|uniq|awk -F \"|\" '\$3/\$2>=0.8 {print \$1}' `; do grep -P \"\$i\" $opt_d/Summary_table.scored.recalled.txt|head -1|awk '\$4>=10 {print}'>>$opt_d/Top_candidate_orthologs.txt; done;
sort -r -nk4 $opt_d/Top_candidate_orthologs.txt >$opt_d/tmp.txt; mv $opt_d/tmp.txt $opt_d/Top_candidate_orthologs.txt;");

}




