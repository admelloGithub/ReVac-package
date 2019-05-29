#!/usr/bin/perl -w

use strict;

my $input=$ARGV[0];
my $cut_off=$ARGV[1];
my $odir=$ARGV[2];
my @inputs;



unless (-e "$odir/SPAAN2table.out") {
`touch $odir/SPAAN2table.out`;
}
open (my $fh3, "<", $input) || die "Couldn't open input list:$!\n";

@inputs = <$fh3>;

close $fh3;

foreach (@inputs) {
    my $file=$_;
    $file =~ s/^\s+//;	
    $file =~ s/\s+$//;

open (my $fh, "<", $file) || die "Couldn't open input:$!\n";

my (@lines,@seq_ids);
@lines = <$fh>;
shift @lines;

close $fh;

foreach (@lines) {

	my @split=split("\t",$_);
	$split[2] =~ s/^\s+//;
	$split[2] =~ s/\s+$//;
	$split[2] =~ s/>//;
	$split[2] =~ s/\.polypeptide\.[0-9]+\.1//;
	if ($split[1] >= $cut_off) {
		push (@seq_ids, "$split[2]\t$split[1]") ;
	}

}

open (my $fh2, ">>", "$odir/SPAAN2table.out") || die "Couldn't open cut off list:$!\n";

foreach (@seq_ids) {
	print $fh2 "$_\n";
}

close $fh2;

}
