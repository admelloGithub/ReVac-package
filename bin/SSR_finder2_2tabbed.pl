#!/usr/bin/perl -w
use strict;
use File::Basename;

my @files;
if ($ARGV[0] =~ /\.list/) {
	my $list = $ARGV[0];
	open (my $fh, "<", $list) || die "$!\n";
        @files = <$fh>;
	close $fh || die "$!\n";
} else {
	@files = @ARGV;
}

foreach (@files) {
    $_ =~ s/\s+$//;
	my ($outlines,$lines);
    my $name = basename($_);
	my $dirname = dirname($_);
	$name =~ s/\..*$//;
	my $fp = $_;
	#$fp =~ s/fsa/fasta/ if ($fp =~ /fsa/);
	$fp =~ s/fsa/fna/ if ($fp =~ /fsa/);
    open (my $fh2, "<", $fp) || die "File:$fp\t$!\n";
		while (<$fh2>) {
			$lines .= $_;
		}
	close $fh2 || die "$!\n";
	
	my @motifs = split(/=/,$lines);
	pop @motifs;

	foreach (@motifs) {
		my @split = split (/\n/,$_);
		my ($pos,$motif,$rep,$seq);		
			foreach my $nuc (@split) {
				if ($nuc =~ /^motif:/) {
					$motif =$nuc;
					$motif =~ s/motif: //;
				} elsif ($nuc =~ /position:/) {
					$pos =$nuc;
					$pos =~ s/position: //;
				} elsif ($nuc =~ /n.rep:/) {
                    $rep =$nuc;
                    $rep =~ s/n.rep: //;
				} elsif ($nuc =~ /seq :/) {
                    $seq =$nuc;
                    $seq =~ s/seq ://;
				}
			}
		$outlines .= "$motif\_$rep\_$seq\t$pos\t$pos\n";
	}
	
	open (my $fh3, ">", "$dirname/$name.ssr") || die "$!\n";
		print $fh3 $outlines if ($lines =~ /=/);
	close $fh3 ||die "$!\n";
}	
