#!/usr/bin/perl -w
use strict;
use List::MoreUtils qw( uniq );
use List::Util qw( min max sum );

#Sample INPUT FILE is ranked table from test_algorithm.pl

my $table = $ARGV[0];
my %seqs;

open (my $fh, "<", $table) || die "$!\n";

##Checking for same number of colums b/w headers and data;
my @headers = split(/\t/,<$fh>);
my @lines = <$fh>;
foreach (@lines) {
	my @linex = split(/\t/,$_);
	if ($_ eq "\t\n") {next;}
	if (scalar(@headers) != scalar(@linex)) {
		die "Improper No of columns between headers and data.\n".scalar(@headers)."\t".scalar(@linex)."\t;$_;";
	}
}
shift @headers;

close $fh || die "$!\n";

my (%cogs,%score,@unmapped);

foreach my $line (@lines) {
    my @split = split(/\t/,$line);
        for (my $i=0; $i<scalar(@headers); $i++) {
			if ($headers[$i] =~ /Jaccard_orthologs/ || $headers[$i] =~ /PanOCT/ || $headers[$i] =~ /Orthomcl/ || $headers[$i] =~ /LSBSR/ || $headers[$i] =~ /Revac_Score/  || $headers[$i] =~ /SSR_Finder/ || $headers[$i] =~ /Attributor/|| $headers[$i] =~ /Length/) {
				$seqs{$split[0]}{$headers[$i]} = $split[$i+1];
			}
		}
}

my @nlines;
foreach my $id (keys %seqs) {
	my $score = $seqs{$id}{Revac_Score};
	$score =~ s/^.*\|//;
	my $nline = "$id\t$score\t";
		foreach my $h (sort keys %{$seqs{$id}}) {
			$nline .= "$seqs{$id}{$h}\t" if ( $h =~ /Jaccard_orthologs/ || $h =~ /PanOCT/ || $h =~ /Orthomcl/ || $h =~ /LSBSR/);
		}
	$nline =~ s/\t$//;
	push (@nlines, $nline);
}

#########################
foreach my $line (@nlines) {

	my @split = split(/\t/,$line);
	$score{$split[0]} = $split[1];
	my (%hash,@keys,$key);
	foreach (@split[2..5]) {
		my $cid = $_;
		$cid =~ s/^.*\|//;
		if ($cid =~ /None/) { next; }
		if ($hash{$cid}) {
			$hash{$cid} .= "$_;#";
    	} else {
			$hash{$cid} = "$_;";
		}	
	}

	foreach (sort keys %hash) {
		if ($hash{$_} =~ /#/) {		#Not trusting single methods with no corroboration from others.
			my $rank;
			if ($hash{$_} =~ /PanOCT/) { $rank += 2; }
			if ($hash{$_} =~ /j_ortholog/) { $rank += 1; }
			if ($hash{$_} =~ /orthomcl/) { $rank += -1; }	#LSBSR is 0
			push (@keys, "$rank\t$hash{$_}"); 
		}
	}
	my $rank;
	foreach (sort @keys) {					#Breaks ties in 2 methods calling different sizes. 
		if ($key) {
			my @split = split(/\t/,$_);
			if ($split[0] > $rank) {
				$key = $split[1];
			} #elsif ($split[0] == $rank) {
				#if ($split[1] =~ /PanOCT/) {  warn "$split[1]\n";$key = $split[1]; }
			#}
		} else {
			my @split = split(/\t/,$_);
			$rank = $split[0];
			$key = $split[1];
		}
	}

	#foreach my $key (@keys) {
		if ($key) {
			$cogs{$key} .= "$split[0],";
		} else {
			push (@unmapped, $line);
		}
	#}
	next;
}

my %ssr;	
my $ct=0;
my @score_filter;
foreach (sort keys %cogs) {
	$ct++;
	my @ct = split(/,/,$cogs{$_});
	my $cnt = scalar(@ct) ;
	my $key = $_;
	$key =~ s/#//g;
	foreach (@ct) {
		push (@score_filter,"$_\t$score{$_}\tOrtholog_$ct|$cnt\t$key\n");
		#print "$_\t$score{$_}\tOrtholog_$ct|$cnt\t$key\n";	#key is all four method ortholog IDs
		$ssr{"Ortholog_$ct|$cnt"} .= "$_;";
	}
}

#print join ("\n",@unmapped); # . scalar(@unmapped) . "\n";


foreach my $oid (sort keys %ssr) {
	$ssr{$oid} =~ s/;$//;
	my @ids = split(/;/,$ssr{$oid});
	my @ssr;
		foreach my $id (@ids) {
			push (@ssr, $seqs{$id}{SSR_Finder}) if ($seqs{$id}{SSR_Finder} =~ /\//);
		}
	my $ssrct = scalar(@ssr);	
	if (scalar(@ids) == scalar(@ssr)) {
		my @ssrs = uniq @ssr;
		#print "$oid\t$ssrct\t@ssrs\n";  #Ortholog_ID,proteins with ssrs, and uniq ssrs.
	}
}

my @filter;
foreach (@score_filter) {
	my @split = split(/\t/,$_);
    $split[1] =~ s/^.*\|//;
    my @scores = split(/-/,$split[1]);
    $scores[0] =~ s/\n$//;
    if ($scores[0] == 0) {
        push (@filter,$_) if ($_ =~ /Ortholog/);
    }
    if ($scores[1]/$scores[0] <= 0.1) {
        push (@filter,$_) if ($_ =~ /Ortholog/);
    }
}

#print @filter;
my %filtered;
foreach (@score_filter) {
	my @split = split(/\t/,$_);
	push(@{$filtered{$split[2]}}, eval($split[1]));
} 
my %ct;
foreach (@filter) {
	my @split = split(/\t/,$_);
	$ct{$split[2]} += 1;
}

print "#LocusTag\tScoreBreakdown\tOrthologID|AllCount|<10perCt\tAllAvg\tAllMin\tAllMax\tAnnotation\tSSR\n";
foreach my $id (@filter) {
	my @split = split(/\t/,$id);
	my $k = $split[2];
	print join ("\t",@split[0..1]) ."\t";
	print "$k|$ct{$k}\t".sum(@{$filtered{$k}})/scalar(@{$filtered{$k}})."\t".min(@{$filtered{$k}})."\t".max(@{$filtered{$k}})."\t".$seqs{$split[0]}{Attributor}."\t".$seqs{$split[0]}{SSR_Finder}."\n";
}



#####
