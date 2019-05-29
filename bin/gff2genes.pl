#!/usr/bin/perl -w
use strict;
use File::Basename;

my $list = $ARGV[0];
my $outdir = $ARGV[1];

open (my $fh, "<", $list) || die "$!\n";
	my @files = <$fh>;
close $fh || die "$!\n";

foreach (@files) {
	$_ =~ s/\s+$//;
	my $name = basename($_);
	$name =~ s/\.gff//;
	#print "$name\n";
	if ($name =~ /_gb[f]+/) { $name =~ s/_gb[f]+//g; }
	my ($contig,%genes);
	open (my $fh2, "<", $_) || die "File:$_\t$!\n";
		while (<$fh2>) {
			if ($_ =~ /=/) { $_ =~ s/=/ /g; }
			if ($_ =~ /;/) { $_ =~ s/;/ ; /g; }
			if ($_ =~ /^#/) { next; }
			elsif ($_ !~ /^#/) { #($_ =~ /^N$name/ || $_ =~ /$name/) {
				my @split = split (/\s+/, $_);
				if ($split[2] eq "CDS" || $split[2] eq "misc_feature") {
					my $lt;
					for (my $i=0; $i<scalar(@split); $i++) {
						if ($split[$i] eq "locus_tag" && $split[6] eq "+") {
							$lt = $split[$i+1];
							$genes{$split[0]}{$split[$i+1]} .= "$split[$i+1]\t$split[3]\t$split[4]\t";
						} elsif ($split[$i] eq "locus_tag" && $split[6] eq "-") {
                            $lt = $split[$i+1];
							$genes{$split[0]}{$split[$i+1]} .= "$split[$i+1]\t$split[4]\t$split[3]\t";
						} elsif ($split[$i] eq "product") {
							for (my $j= $i+1; ($split[$j] ne ";" || $split[$j] =~ /\"$/); $j++) {
								$genes{$split[0]}{$lt} .= "$split[$j] ";
							}
						}
					}
				} elsif ($split[2] eq "tRNA") {
					next;
					my $lt;
                    for (my $i=8; $i<scalar(@split); $i++) {
                        if ($split[$i] eq "locus_tag" && $split[6] eq "+") {
                            $lt = $split[$i+1];
                            $genes{$split[0]}{$split[$i+1]} .= "$split[$i+1]\t$split[3]\t$split[4]\t";
                        } elsif ($split[$i] eq "locus_tag" && $split[6] eq "-") {
                            $lt = $split[$i+1];
                            $genes{$split[0]}{$split[$i+1]} .= "$split[$i+1]\t$split[4]\t$split[3]\t";
                        } elsif ($split[$i] eq "product" && $split[$i+1] =~ /^"tRNA/) {
                            for (my $j= $i+1; $split[$j] ne ";"; $j++) {
                                $genes{$split[0]}{$lt} .= "$split[$j] ";
                            }
                        }
                    }
				}
			}
		}
	close $fh2 || die "$!\n";

	foreach my $contig (sort keys %genes) {
		my $file = $name;
		$file =~ s/\..*$//;
		open (my $fh4, ">", "$outdir/$file.genes") || die "$!\n";
		foreach (sort keys %{$genes{$contig}} ) {
			my $value = $genes{$contig}{$_};
			$value =~ s/"//g;
			print $fh4 "$value\n";
		}
		close $fh4 || die "$!\n";
	}

}
