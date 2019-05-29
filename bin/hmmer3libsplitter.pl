#!usr/bin/perl -w

use strict;

my $hmm_file = $ARGV[0];
my $desc = $ARGV[1];
my $out_file = $ARGV[2];
my @lines;
my @proteins;
my @pls;

open(my $fh,"<", $hmm_file) or die "Couldn't open hmmfile:$!\n";

	@lines = <$fh>;

close $fh;

my $count=0;

foreach(@lines) {

	$count++;	
	if ( $_ =~/$desc/ ) {
		
	my $i=$count-4;
		until($lines[$i] =~ /\/\//) {
		
			push @pls, $lines[$i];
			$i++;
		}
	push @pls, "//\n";
	}
}

open(my $fh2, ">", $out_file) or die "Couldn't open outputfile:$!\n";

	foreach(@pls) {
		#$_ =~ s/^\s+//;
		print $fh2 "$_";
	}

close $fh2;
