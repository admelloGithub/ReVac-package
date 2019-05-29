#! /usr/bin/perl -w

use File::Basename;
#use lib '/home/sienaem1/new_ssr_search/';
#use  Text::LevenshteinXS qw(distance);
use lib '/local/projects/PNTHI/herve_dir/ssr_dir/SSR_finder/';
use Text::WagnerFischer qw(distance);


#####Originally designed for DNA tandem repeats, modified to locate protein repeat sequences. Allows mismatches.


####USAGE
#1st argument is an output directory followed by a list of files or just space seperated file paths.
# eg:  ./SSR_finder_genes.pl ./ file1 file2 .... 

# get the list of fasta files
#@ARGV = glob ("*.{fna,fasta}");
my $outdir = $ARGV[0];
shift @ARGV;

if ($ARGV[0] =~ /\.list/) {
    my $input_list = $ARGV[1];
    open (my $input, "<", $input_list) || die "Couldn't open input individual files list: $!\n";
    @ARGV = <$input>;
    close $input || die "$!\n";
}

# loop over each file and process it
foreach $dnasequence(@ARGV) {

	chomp $dnasequence;
	
	#check if the specified file exists
	unless (-e $dnasequence){
	       print "File \"$dnasequence\" doesn't seem to exist!\n\n";
	       exit;   
	}
	
	#check if it's possible to open the specified file
	unless (open(DNAFILE, $dnasequence)){
	       print "Cannot open file \"$dnasequence\"\n\n";
	       exit;
	}
	
	#store the fileinto an array
	my @DNA = <DNAFILE>;
	close DNAFILE;
	my ($id,%seqs);
	#remove the first line of text (the first line in a file fasta is the heading)
	foreach (@DNA){
		if ($_ =~ />/ ){
		$id = $_;
		} else {
	    $seqs{$id} .= $_;
	    }
	}

	($outputfile) = basename($dnasequence); #=~ /^(\w+)\.\w+$/;
    $outputfile = "$outdir/$outputfile.ssr.txt";

	#open a filehandle for the output
    unless (-e $outputfile) {
    open STDOUT, ">> $outputfile";
    } else {
    system ("rm $outputfile");
    open STDOUT, ">> $outputfile";
    }

	
	foreach (keys %seqs) {
	#save the sequence in a single string and remove empty spaces and newlines
	#$DNA = join ("",@DNA);
	$DNA = $seqs{$_};
	$DNA =~ s/[ \n]//g;
	$DNA = uc $DNA;
	
	#specify the length of the DNA sequence to be analyzed
	$numero_delle_basi = length $DNA;

	#ask the name for an output file and open the relative filehandle
	# print "\nPlease type the name for the output file: ";
	#($outputfile) = $dnasequence =~ /^(\w+)\.\w+$/;
	#($outputfile) = basename($dnasequence); #=~ /^(\w+)\.\w+$/;
	#$outputfile = "$outdir/$outputfile.ssr.txt";
	
	#print some preliminary information at the beginning of the output file
	print "\n\nfile analyzed: $dnasequence\n\n";
	print "$_\n";
	print "Sequence length: $numero_delle_basi\n\n";
	
	#initialize the counts of the motifs
	my $number_tandem_repeats_1 = 0;
	my $number_tandem_repeats_2 = 0;
	my $number_tandem_repeats_3 = 0;
	my $number_tandem_repeats_4 = 0;
	my $number_tandem_repeats_5 = 0;
	my $number_tandem_repeats_6 = 0;
	my $number_tandem_repeats_7 = 0;
	my $number_tandem_repeats_8 = 0;
	my $number_tandem_repeats_9 = 0;
	my $number_tandem_repeats_10= 0;

	#this is a loop that will scan the DNA sequence for repeated motifs. ten different conditional sub-loops set different parameters depending on the length of the motif analyzed (from 1 to 10 bases).
	for ( my $Position = 0 ; $Position < length $DNA ; $Position=($Position + 1 + $jump)) {
		$st = substr ($DNA, $Position, 200);
		$jump=0;

	#>>>>>> 1
	$motif = substr ($st, 0, 1);                  #change here!
		@elements = &splice_($st, 1);         #change here!
		$t = 0;
		$number_repeats = 0;
		foreach $r(@elements){
			if (distance([0,5,1], $motif, $r) < 1){  #change here!
				$number_repeats += 1;
			}elsif (distance([0,5,1], $motif, $r) < 1 and $t < 1 ){  #change here!
				$number_repeats += 1;
				$t++;
			}else{
				last;
			}
		}
		
		undef @elements;
		$SSR = substr ($st, 0, (1 * $number_repeats));  #change here!

		
		if ($number_repeats > 8 and $SSR =~ /(\w{1})\1{2,}/){           #decide the minimum number of repeats #HT Moxon
			if ($SSR !~ /$motif$/){
				substr ($SSR, -1) = "";                #change here!
				$number_repeats--;
			}
			print "motif: $motif\nposition: ", $Position +1 ,"\nn.rep: $number_repeats\nseq :$SSR\n=\n";
			$jump = ((length $motif) * $number_repeats) - 1;
			$number_tandem_repeats_1++;
			next;
		}
	#<<<<<< 1


	#>>>>>> 2
	$motif = substr ($st, 0, 2);                  #change here!
		@elements = &splice_($st, 2);         #change here!
		$t = 0;
		$number_repeats = 0;
		foreach $r(@elements){
			if (distance([0,5,1], $motif, $r) < 1){  #change here!
				$number_repeats += 1;
			}elsif (distance([0,5,1], $motif, $r) < 1 and $t < 1 ){  #change here!
				$number_repeats += 1;
				$t++;
			}else{
				last;
			}
		}
		
		undef @elements;
		$SSR = substr ($st, 0, (2 * $number_repeats));  #change here!
		if ($SSR !~ /$motif$/){
			substr ($SSR, -2) = "";                #change here!
			$number_repeats--;
		}
		
		if ($number_repeats > 3 and $SSR !~ /(^\w{1})\1/){            #decide the minimum number of repeats
			print "motif: $motif\nposition: ", $Position +1 ,"\nn.rep: $number_repeats\nseq :$SSR\n=\n";
			$jump = ((length $motif) * $number_repeats) - 1;
			$number_tandem_repeats_2++;
			next;
		}
	#<<<<<< 2


	#>>>>>> 3
	$motif = substr ($st, 0, 3);                  #change here!
		@elements = &splice_($st, 3);         #change here!
		$t = 0;
		$number_repeats = 0;
		foreach $r(@elements){
			if (distance([0,5,1], $motif, $r) < 1){  #change here!
				$number_repeats += 1;
			}elsif (distance([0,5,1], $motif, $r) < 1 and $t < 1 ){  #change here!
				$number_repeats += 1;
				$t++;
			}else{
				last;
			}
		}
		
		undef @elements;
		$SSR = substr ($st, 0, (3 * $number_repeats));  #change here!
		if ($SSR !~ /$motif$/){
			substr ($SSR, -3) = "";                #change here!
			$number_repeats--;
		}
		
		if ($number_repeats > 2 and $SSR !~ /(^\w{1})\1\1/){            #decide the minimum number of repeats
			print "motif: $motif\nposition: ", $Position +1 ,"\nn.rep: $number_repeats\nseq :$SSR\n=\n";
			$jump = ((length $motif) * $number_repeats) - 1;
			$number_tandem_repeats_3++;
			next;
		}
	#<<<<<< 3


	#>>>>>> 4
	$motif = substr ($st, 0, 4);                  #change here!
		@elements = &splice_($st, 4);         #change here!
		$t = 0;
		$number_repeats = 0;
		foreach $r(@elements){
			if (distance([0,5,1], $motif, $r) < 1){  #change here!
				$number_repeats += 1;
			}elsif (distance([0,5,1], $motif, $r) < 1 and $t < 1 ){  #change here!
				$number_repeats += 1;
				$t++;
			}else{
				last;
			}
		}
		
		undef @elements;
		$SSR = substr ($st, 0, (4 * $number_repeats));  #change here!
	# 	if ($SSR !~ /$motif$/){
	# 		substr ($SSR, -4) = "";                #change here!
	# 		$number_repeats--;
	# 	}
		
		if ($number_repeats > 2 and $SSR !~ /(^\w{1})\1\1\1/) {  #and $SSR !~ /(TTAA){3}/){            #decide the minimum number of repeats
			print "motif: $motif\nposition: ", $Position +1 ,"\nn.rep: $number_repeats\nseq :$SSR\n=\n";
			$jump = ((length $motif) * $number_repeats) - 1;
			$number_tandem_repeats_4++;
			next;
		}
	#<<<<<< 4


	#>>>>>> 5
	$motif = substr ($st, 0, 5);                  #change here!
		@elements = &splice_($st, 5);         #change here!
		$t = 0;
		$number_repeats = 0;
		foreach $r(@elements){
			if (distance([0,5,1], $motif, $r) < 1){  #change here!
				$number_repeats += 1;
			}elsif (distance([0,5,1], $motif, $r) < 2 and $t < 1 ){  #change here!
				$number_repeats += 1;
				$t++;
			}else{
				last;
			}
		}
		
		undef @elements;
		$SSR = substr ($st, 0, (5 * $number_repeats));  #change here!
	# 	if ($SSR !~ /$motif$/){
	# 		substr ($SSR, -5) = "";                #change here!
	# 		$number_repeats--;
	# 	}
		
# 		if ($number_repeats > 3)            #decide the minimum number of repeats #HT Moxon
		if ($number_repeats > 2){            #decide the minimum number of repeats #HT Moxon
			print "motif: $motif\nposition: ", $Position +1 ,"\nn.rep: $number_repeats\nseq :$SSR\n=\n";
			$jump = ((length $motif) * $number_repeats) - 1;
			$number_tandem_repeats_5++;
			next;
		}
	#<<<<<< 5


	#>>>>>> 6
	$motif = substr ($st, 0, 6);                  #change here!
		@elements = &splice_($st, 6);         #change here!
		$t = 0;
		$number_repeats = 0;
		foreach $r(@elements){
			if (distance([0,5,1], $motif, $r) < 1){  #change here!
				$number_repeats += 1;
			}elsif (distance([0,5,1], $motif, $r) < 2 and $t < 1 ){  #change here!
				$number_repeats += 1;
				$t++;
			}else{
				last;
			}
		}
		
		undef @elements;
		$SSR = substr ($st, 0, (6 * $number_repeats));  #change here!
	# 	if ($SSR !~ /$motif$/){
	# 		substr ($SSR, -6) = "";                #change here!
	# 		$number_repeats--;
	# 	}
		
		if ($number_repeats >= 2){            #decide the minimum number of repeats
			print "motif: $motif\nposition: ", $Position +1 ,"\nn.rep: $number_repeats\nseq :$SSR\n=\n";
			$jump = ((length $motif) * $number_repeats) - 1;
			$number_tandem_repeats_6++;
			next;
		}
	#<<<<<< 6


	#>>>>>> 7
	$motif = substr ($st, 0, 7);                  #change here!
		@elements = &splice_($st, 7);         #change here!
		$t = 0;
		$number_repeats = 0;
		foreach $r(@elements){
			if (distance([0,5,1], $motif, $r) < 1){  #change here!
				$number_repeats += 1;
			}elsif (distance([0,5,1], $motif, $r) < 2 and $t < 1 ){  #change here!
				$number_repeats += 1;
				$t++;
			}else{
				last;
			}
		}
		
		undef @elements;
		$SSR = substr ($st, 0, (7 * $number_repeats));  #change here!
	# 	if ($SSR !~ /$motif$/){
	# 		substr ($SSR, -7) = "";                #change here!
	# 		$number_repeats--;
	# 	}
		
		if ($number_repeats >= 2){            #decide the minimum number of repeats
			print "motif: $motif\nposition: ", $Position +1 ,"\nn.rep: $number_repeats\nseq :$SSR\n=\n";
			$jump = ((length $motif) * $number_repeats) - 1;
			$number_tandem_repeats_7++;
			next;
		}
	#<<<<<< 7


	#>>>>>> 8
	$motif = substr ($st, 0, 8);                  #change here!
		@elements = &splice_($st, 8);         #change here!
		$t = 0;
		$number_repeats = 0;
		foreach $r(@elements){
			if (distance([0,5,1], $motif, $r) < 1){  #change here!
				$number_repeats += 1;
			}elsif (distance([0,5,1], $motif, $r) < 3 and $t < 1 ){  #change here!
				$number_repeats += 1;
				$t++;
			}else{
				last;
			}
		}
		
		undef @elements;
		$SSR = substr ($st, 0, (8 * $number_repeats));  #change here!
	# 	if ($SSR !~ /$motif$/){
	# 		substr ($SSR, -8) = "";                #change here!
	# 		$number_repeats--;
	# 	}
		
		if ($number_repeats >= 2){            #decide the minimum number of repeats
			print "motif: $motif\nposition: ", $Position +1 ,"\nn.rep: $number_repeats\nseq :$SSR\n=\n";
			$jump = ((length $motif) * $number_repeats) - 1;
			$number_tandem_repeats_8++;
			next;
		}
	#<<<<<< 8



	#>>>>>> 9
	$motif = substr ($st, 0, 9);                  #change here!
		@elements = &splice_($st, 9);         #change here!
		$t = 0;
		$number_repeats = 0;
		foreach $r(@elements){
			if (distance([0,5,1], $motif, $r) < 1){  #change here!
				$number_repeats += 1;
			}elsif (distance([0,5,1], $motif, $r) < 3 and $t < 1 ){  #change here!
				$number_repeats += 1;
				$t++;
			}else{
				last;
			}
		}
		
		undef @elements;
		$SSR = substr ($st, 0, (9 * $number_repeats));  #change here!
	# 	if ($SSR !~ /$motif$/){
	# 		substr ($SSR, -9) = "";                #change here!
	# 		$number_repeats--;
	# 	}
		
		if ($number_repeats >= 2){            #decide the minimum number of repeats
			print "motif: $motif\nposition: ", $Position +1 ,"\nn.rep: $number_repeats\nseq :$SSR\n=\n";
			$jump = ((length $motif) * $number_repeats) - 1;
			$number_tandem_repeats_9++;
			next;
		}
	#<<<<<< 9


	#>>>>>> 10
	$motif = substr ($st, 0, 10);                  #change here!
		@elements = &splice_($st, 10);         #change here!
		$t = 0;
		$number_repeats = 0;
		foreach $r(@elements){
			if (distance([0,5,1], $motif, $r) < 1){  #change here!
				$number_repeats += 1;
			}elsif (distance([0,5,1], $motif, $r) < 3 and $t < 1 ){  #change here!
				$number_repeats += 1;
				$t++;
			}else{
				last;
			}
		}
		
		undef @elements;
		$SSR = substr ($st, 0, (10 * $number_repeats));  #change here!
	# 	if ($SSR !~ /$motif$/){
	# 		substr ($SSR, -10) = "";                #change here!
	# 		$number_repeats--;
	# 	}
		
		if ($number_repeats >= 2){            #decide the minimum number of repeats
			print "motif: $motif\nposition: ", $Position +1 ,"\nn.rep: $number_repeats\nseq :$SSR\n=\n";
			$jump = ((length $motif) * $number_repeats) - 1;
			$number_tandem_repeats_10++;
			next;
		}
	#<<<<<< 10

	#>>>>>> 11
     $motif = substr ($st, 0, 11);                  #change here!
         @elements = &splice_($st, 11);         #change here!
         $t = 0;
         $number_repeats = 0;
         foreach $r(@elements){
             if (distance([0,5,1], $motif, $r) < 1){  #change here!
                 $number_repeats += 1;
             }elsif (distance([0,5,1], $motif, $r) < 3 and $t < 1 ){  #change here!
                 $number_repeats += 1;
                 $t++;
             }else{
                 last;
             }
         }
 
         undef @elements;
         $SSR = substr ($st, 0, (11 * $number_repeats));  #change here!
     #   if ($SSR !~ /$motif$/){
     #       substr ($SSR, -10) = "";                #change here!
     #       $number_repeats--;
     #   }
 
         if ($number_repeats >= 2){            #decide the minimum number of repeats
             print "motif: $motif\nposition: ", $Position +1 ,"\nn.rep: $number_repeats\nseq :$SSR\n=\n";
             $jump = ((length $motif) * $number_repeats) - 1;
             $number_tandem_repeats_10++;
             next;
         }
     #<<<<<< 11

	#>>>>>> 12
     $motif = substr ($st, 0, 12);                  #change here!
         @elements = &splice_($st, 12);         #change here!
         $t = 0;
         $number_repeats = 0;
         foreach $r(@elements){
             if (distance([0,5,1], $motif, $r) < 1){  #change here!
                 $number_repeats += 1;
             }elsif (distance([0,5,1], $motif, $r) < 3 and $t < 1 ){  #change here!
                 $number_repeats += 1;
                 $t++;
             }else{
                 last;
             }
         }
 
         undef @elements;
         $SSR = substr ($st, 0, (12 * $number_repeats));  #change here!
     #   if ($SSR !~ /$motif$/){
     #       substr ($SSR, -10) = "";                #change here!
     #       $number_repeats--;
     #   }
 
         if ($number_repeats >= 2){            #decide the minimum number of repeats
             print "motif: $motif\nposition: ", $Position +1 ,"\nn.rep: $number_repeats\nseq :$SSR\n=\n";
             $jump = ((length $motif) * $number_repeats) - 1;
             $number_tandem_repeats_10++;
             next;
         }
     #<<<<<< 12	

	#>>>>>> 13
     $motif = substr ($st, 0, 13);                  #change here!
         @elements = &splice_($st, 13);         #change here!
         $t = 0;
         $number_repeats = 0;
         foreach $r(@elements){
             if (distance([0,5,1], $motif, $r) < 1){  #change here!
                 $number_repeats += 1;
             }elsif (distance([0,5,1], $motif, $r) < 3 and $t < 1 ){  #change here!
                 $number_repeats += 1;
                 $t++;
             }else{
                 last;
             }
         }
 
         undef @elements;
         $SSR = substr ($st, 0, (13 * $number_repeats));  #change here!
     #   if ($SSR !~ /$motif$/){
     #       substr ($SSR, -10) = "";                #change here!
     #       $number_repeats--;
     #   }
 
         if ($number_repeats >= 2){            #decide the minimum number of repeats
             print "motif: $motif\nposition: ", $Position +1 ,"\nn.rep: $number_repeats\nseq :$SSR\n=\n";
             $jump = ((length $motif) * $number_repeats) - 1;
             $number_tandem_repeats_10++;
             next;
         }
     #<<<<<< 13

	#>>>>>> 14
     $motif = substr ($st, 0, 14);                  #change here!
         @elements = &splice_($st, 14);         #change here!
         $t = 0;
         $number_repeats = 0;
         foreach $r(@elements){
             if (distance([0,5,1], $motif, $r) < 1){  #change here!
                 $number_repeats += 1;
             }elsif (distance([0,5,1], $motif, $r) < 3 and $t < 1 ){  #change here!
                 $number_repeats += 1;
                 $t++;
             }else{
                 last;
             }
         }
 
         undef @elements;
         $SSR = substr ($st, 0, (14 * $number_repeats));  #change here!
     #   if ($SSR !~ /$motif$/){
     #       substr ($SSR, -10) = "";                #change here!
     #       $number_repeats--;
     #   }
 
         if ($number_repeats >= 2){            #decide the minimum number of repeats
             print "motif: $motif\nposition: ", $Position +1 ,"\nn.rep: $number_repeats\nseq :$SSR\n=\n";
             $jump = ((length $motif) * $number_repeats) - 1;
             $number_tandem_repeats_10++;
             next;
         }
     #<<<<<< 14

	#>>>>>> 15
     $motif = substr ($st, 0, 15);                  #change here!
         @elements = &splice_($st, 15);         #change here!
         $t = 0;
         $number_repeats = 0;
         foreach $r(@elements){
             if (distance([0,5,1], $motif, $r) < 1){  #change here!
                 $number_repeats += 1;
             }elsif (distance([0,5,1], $motif, $r) < 3 and $t < 1 ){  #change here!
                 $number_repeats += 1;
                 $t++;
             }else{
                 last;
             }
         }
 
         undef @elements;
         $SSR = substr ($st, 0, (15 * $number_repeats));  #change here!
     #   if ($SSR !~ /$motif$/){
     #       substr ($SSR, -10) = "";                #change here!
     #       $number_repeats--;
     #   }
 
         if ($number_repeats >= 2){            #decide the minimum number of repeats
             print "motif: $motif\nposition: ", $Position +1 ,"\nn.rep: $number_repeats\nseq :$SSR\n=\n";
             $jump = ((length $motif) * $number_repeats) - 1;
             $number_tandem_repeats_10++;
             next;
         }
     #<<<<<< 15

	#>>>>>> 16
     $motif = substr ($st, 0, 16);                  #change here!
         @elements = &splice_($st, 16);         #change here!
         $t = 0;
         $number_repeats = 0;
         foreach $r(@elements){
             if (distance([0,5,1], $motif, $r) < 1){  #change here!
                 $number_repeats += 1;
             }elsif (distance([0,5,1], $motif, $r) < 3 and $t < 1 ){  #change here!
                 $number_repeats += 1;
                 $t++;
             }else{
                 last;
             }
         }
 
         undef @elements;
         $SSR = substr ($st, 0, (16 * $number_repeats));  #change here!
     #   if ($SSR !~ /$motif$/){
     #       substr ($SSR, -10) = "";                #change here!
     #       $number_repeats--;
     #   }
 
         if ($number_repeats >= 2){            #decide the minimum number of repeats
             print "motif: $motif\nposition: ", $Position +1 ,"\nn.rep: $number_repeats\nseq :$SSR\n=\n";
             $jump = ((length $motif) * $number_repeats) - 1;
             $number_tandem_repeats_10++;
             next;
         }
     #<<<<<< 16

	#>>>>>> 17
     $motif = substr ($st, 0, 17);                  #change here!
         @elements = &splice_($st, 17);         #change here!
         $t = 0;
         $number_repeats = 0;
         foreach $r(@elements){
             if (distance([0,5,1], $motif, $r) < 1){  #change here!
                 $number_repeats += 1;
             }elsif (distance([0,5,1], $motif, $r) < 3 and $t < 1 ){  #change here!
                 $number_repeats += 1;
                 $t++;
             }else{
                 last;
             }
         }
 
         undef @elements;
         $SSR = substr ($st, 0, (17 * $number_repeats));  #change here!
     #   if ($SSR !~ /$motif$/){
     #       substr ($SSR, -10) = "";                #change here!
     #       $number_repeats--;
     #   }
 
         if ($number_repeats >= 2){            #decide the minimum number of repeats
             print "motif: $motif\nposition: ", $Position +1 ,"\nn.rep: $number_repeats\nseq :$SSR\n=\n";
             $jump = ((length $motif) * $number_repeats) - 1;
             $number_tandem_repeats_10++;
             next;
         }
     #<<<<<< 17

	#>>>>>> 18
     $motif = substr ($st, 0, 18);                  #change here!
         @elements = &splice_($st, 18);         #change here!
         $t = 0;
         $number_repeats = 0;
         foreach $r(@elements){
             if (distance([0,5,1], $motif, $r) < 1){  #change here!
                 $number_repeats += 1;
             }elsif (distance([0,5,1], $motif, $r) < 3 and $t < 1 ){  #change here!
                 $number_repeats += 1;
                 $t++;
             }else{
                 last;
             }
         }
 
         undef @elements;
         $SSR = substr ($st, 0, (18 * $number_repeats));  #change here!
     #   if ($SSR !~ /$motif$/){
     #       substr ($SSR, -10) = "";                #change here!
     #       $number_repeats--;
     #   }
 
         if ($number_repeats >= 2){            #decide the minimum number of repeats
             print "motif: $motif\nposition: ", $Position +1 ,"\nn.rep: $number_repeats\nseq :$SSR\n=\n";
             $jump = ((length $motif) * $number_repeats) - 1;
             $number_tandem_repeats_10++;
             next;
         }
     #<<<<<< 18

	#>>>>>> 19
     $motif = substr ($st, 0, 19);                  #change here!
         @elements = &splice_($st, 19);         #change here!
         $t = 0;
         $number_repeats = 0;
         foreach $r(@elements){
             if (distance([0,5,1], $motif, $r) < 1){  #change here!
                 $number_repeats += 1;
             }elsif (distance([0,5,1], $motif, $r) < 3 and $t < 1 ){  #change here!
                 $number_repeats += 1;
                 $t++;
             }else{
                 last;
             }
         }
 
         undef @elements;
         $SSR = substr ($st, 0, (19 * $number_repeats));  #change here!
     #   if ($SSR !~ /$motif$/){
     #       substr ($SSR, -10) = "";                #change here!
     #       $number_repeats--;
     #   }
 
         if ($number_repeats >= 2){            #decide the minimum number of repeats
             print "motif: $motif\nposition: ", $Position +1 ,"\nn.rep: $number_repeats\nseq :$SSR\n=\n";
             $jump = ((length $motif) * $number_repeats) - 1;
             $number_tandem_repeats_10++;
             next;
         }
     #<<<<<< 19

	#>>>>>> 20
     $motif = substr ($st, 0, 20);                  #change here!
         @elements = &splice_($st, 20);         #change here!
         $t = 0;
         $number_repeats = 0;
         foreach $r(@elements){
             if (distance([0,5,1], $motif, $r) < 1){  #change here!
                 $number_repeats += 1;
             }elsif (distance([0,5,1], $motif, $r) < 3 and $t < 1 ){  #change here!
                 $number_repeats += 1;
                 $t++;
             }else{
                 last;
             }
         }
 
         undef @elements;
         $SSR = substr ($st, 0, (20 * $number_repeats));  #change here!
     #   if ($SSR !~ /$motif$/){
     #       substr ($SSR, -10) = "";                #change here!
     #       $number_repeats--;
     #   }
 
         if ($number_repeats >= 2){            #decide the minimum number of repeats
             print "motif: $motif\nposition: ", $Position +1 ,"\nn.rep: $number_repeats\nseq :$SSR\n=\n";
             $jump = ((length $motif) * $number_repeats) - 1;
             $number_tandem_repeats_10++;
             next;
         }
     #<<<<<< 20

	}



	print "\nIdentified SSRs ( 1b): $number_tandem_repeats_1\n";
	print "Identified SSRs ( 2b): $number_tandem_repeats_2\n";
	print "Identified SSRs ( 3b): $number_tandem_repeats_3\n";
	print "Identified SSRs ( 4b): $number_tandem_repeats_4\n";
	print "Identified SSRs ( 5b): $number_tandem_repeats_5\n";
	print "Identified SSRs ( 6b): $number_tandem_repeats_6\n";
	print "Identified SSRs ( 7b): $number_tandem_repeats_7\n";
	print "Identified SSRs ( 8b): $number_tandem_repeats_8\n";
	print "Identified SSRs ( 9b): $number_tandem_repeats_9\n";
	print "Identified SSRs (10b): $number_tandem_repeats_10\n\n";

	my $number_tandem_repeats_tot =
	   $number_tandem_repeats_1 +
	   $number_tandem_repeats_2 +
	   $number_tandem_repeats_3 +
	   $number_tandem_repeats_4 +
	   $number_tandem_repeats_5 +
	   $number_tandem_repeats_6 +
	   $number_tandem_repeats_7 +
	   $number_tandem_repeats_8 +
	   $number_tandem_repeats_9 +
	   $number_tandem_repeats_10;

	print "Identified SSRs (TOT): $number_tandem_repeats_tot\n\n";

	#close STDOUT;

	}

	close STDOUT;	

	}


#exit;

sub similarity {
	(my $query, my $target, my $treshold) = @_;
	$score = distance ($query, $target);
	if ($score >= $treshold){
		return 1;
	}else{
		return 0;
	}
}

sub splice_ {
	my $string = $_[0];
	my $size = $_[1];
	my @array = ();
	until (length $string == 0){
		push @array , (substr ($string, 0, $size, "") );
	}
	return @array;
}
