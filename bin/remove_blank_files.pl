#!usr/bin/perl -w
use strict;

my $i; 
my @pttlines;
my $k;
my $pl;

foreach (@ARGV) {
	
	if ( $_ =~ /.ptt/ ) {
	my $pttfile = $_;	
	my $pttout = $_ . ".fixed.list";
	system ("touch $pttout");
	open (my $fh3, "<", $pttfile) or die "Couldn't open ptt list file : $!\n";

                @pttlines = <$fh3>;

        close $fh3;
        $pl =scalar(@pttlines);
	
	open (my $fh5, ">", $pttout) or die "Couldn't open ptt fixed list: $!\n";

	for ($k=0; $k < $pl; $k++ ) {	
		#print "$pttlines[$k]\n";	
		$pttlines[$k] =~ s/^\s+//;
		$pttlines[$k] =~ s/\s+$//;

		open (my $fh4, "<", $pttlines[$k]) or die "Couldn't open ptt file : $! $pttlines[$k]\n";

		my $pttline1 = <$fh4>;
		my $pttline2 = <$fh4>;
		
			if (!( $pttline2 =~ /^0 proteins$/ )) {

			print $fh5 "$pttlines[$k]\n";

			}
	
		close $fh4;
	}

	close $fh5;
	system ("mv $pttout $pttfile");

	}

	else {
	my $file = $_;
	my $out = $_ . ".fixed.list";

	my @lines;
	my @newlines;	

	system ("touch $out");
	open (my $fh, "<", $file) or die "Couldn't open list file : $!\n";
	
		@lines = <$fh>;
			
	close $fh;
	my $l =scalar(@lines);
	
	for ( $i=0; $i < $l; $i++ ) {
              
		$lines[$i] =~ s/^\s+//;
		$lines[$i] =~ s/\s+$//;
                if (!(-z $lines[$i])) {
              		push (@newlines, $lines[$i]);
                }

	}

        open (my $fh2, ">", $out) or die "Couldn't open list file : $!\n";

        for ( my $j=0; $j < scalar(@newlines); $j++) {
		
		print $fh2 "$newlines[$j]\n";
	
	}

	close $fh2;

	system ("mv $out $file");

	}
}
