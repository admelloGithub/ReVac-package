#!/usr/local/bin/perl -w
use strict;

foreach (@ARGV) {

my $multifasta = $_;
my %seqs;

open (my $fh, "<", $multifasta) || die "Couldn't open input multifasta file for reading:$!\n";
    my $text = do { local $/; <$fh> };
close $fh || die "$!\n";

my @seqs = split (/>/, $text);
shift @seqs;

foreach (@seqs) {
    my @split = split (/\n/, $_);
    $split[0] =~ s/^/>/;
    $seqs{$split[0]} = join('', @split[1..(scalar(@split)-1)]);
}

if (-e "$multifasta.trimmed") {exit (0);}

open (my $fh2, ">", "$multifasta.trimmed") || die "Couldn't open input multifasta file for writing:$!\n";

foreach (keys %seqs) {
    if ($seqs{$_} =~ /X+/) {
			#print $fh2 "$_\n";
		my @parts = split(/X+/, $seqs{$_});
		if (scalar(@parts)!=2) {
			$seqs{$_} = $parts[0];
		} else {

			if (length($parts[0]) > length($parts[1])) {
				$seqs{$_} = $parts[0];
			}    
			elsif (length($parts[0]) < length($parts[1])) {        
				$seqs{$_} = $parts[1];
			}
			elsif (length($parts[0]) == length ($parts[1])) {
				$seqs{$_} = $parts[0];
			}
		}

    }
	$seqs{$_} =~ s/(.{1,60})/$1\n/gs;
	$seqs{$_} =~ s/(.{1,})/$1/gs;
	print $fh2 "$_\n$seqs{$_}" if (length($seqs{$_}) >15);
}		

close $fh2 || die "$!\n";

}
