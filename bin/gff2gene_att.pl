#!/usr/bin/perl -w
use strict;
use File::Basename;
use List::MoreUtils qw(uniq);

my $list = $ARGV[0];
my $outdir = $ARGV[1];
my %genes;
open (my $fh, "<", $list) || die "$!\n";
	my @files = <$fh>;
close $fh || die "$!\n";

foreach (@files) {
	$_ =~ s/\s+$//;
	my $name = basename($_);
	$name =~ s/\.gff//;
	open (my $fh2, "<", $_) || die "File:$_\t$!\n";
		while (<$fh2>) {
			if ($_ =~ /^#/) { next; }
			else { #if ($_ =~ /^N$name/) {
				my $line = $_;
				$line =~ s/=/ =/g;
				$line =~ s/;/ ; /g;
				my @split = split (/\s+/, $line);
				if ($split[2] eq "CDS" && $line =~ /; translation =/) {  #|| $split[2] eq "misc_feature") {
					my $lt;
					for (my $i=8; $i<scalar(@split); $i++) {
						if ($split[$i] eq "locus_tag" && $split[6] eq "+") {
							$lt = $split[$i+1];
							$lt =~ s/=//g;
							$genes{$split[0]}{$lt} = "$lt\t$split[3]\t$split[4]\t";
						} elsif ($split[$i] eq "locus_tag" && $split[6] eq "-") {
                            $lt = $split[$i+1];
							$lt =~ s/=//g;
							$genes{$split[0]}{$lt} = "$lt\t$split[4]\t$split[3]\t";
						} elsif ($split[$i] eq "product") {
							for (my $j= $i+1; $split[$j] ne ";"; $j++) {
								$genes{$split[0]}{$lt} .= "$split[$j] ";
								$genes{$split[0]}{$lt} =~ s/=//g;
							}
						}
					}
				}
#				} elsif ($split[2] eq "tRNA") {
#					my $lt;
#                    for (my $i=8; $i<scalar(@split); $i++) {
#                        if ($split[$i] eq "locus_tag" && $split[6] eq "+") {
#                            $lt = $split[$i+1];
#                            $genes{$split[0]}{$split[$i+1]} .= "$split[$i+1]\t$split[3]\t$split[4]\t";
#                        } elsif ($split[$i] eq "locus_tag" && $split[6] eq "-") {
#                            $lt = $split[$i+1];
#                            $genes{$split[0]}{$split[$i+1]} .= "$split[$i+1]\t$split[4]\t$split[3]\t";
#                        } elsif ($split[$i] eq "product" && $split[$i+1] =~ /^"tRNA/) {
#                            for (my $j= $i+1; $split[$j] ne ";"; $j++) {
#                                $genes{$split[0]}{$lt} .= "$split[$j] ";
#                            }
#                        }
#                    }
#				}
			}
		}
	close $fh2 || die "$!\n";

}

open (my $fh3, ">", "$outdir/gene_atts.txt") || die "$!\n";
open (my $fh4, ">", "$outdir/gene_tags.txt") || die "$!\n";
my @tags;

foreach my $contig (sort keys %genes) {
	foreach (sort keys %{$genes{$contig}} ) {
		if ($contig =~ /_[0-9]{1,3}\.{0,1}$/) {
			my @split = split(/_/,$contig);
			my $value = $genes{$contig}{$_};
			$value =~ s/"//g;
			#print $fh3 "$split[1]\t$value\t$split[0]\n";
			print $fh3 "$contig\t$value\t$contig\n"; #$plit[0]\n";
			push (@tags, $contig); #$split[0]);
		} else {
			my @split = split(/0{2,}/,$contig);
			my $value = $genes{$contig}{$_};
			$value =~ s/"//g;
			#print $fh3 "$split[$#split]\t$value\t$split[0]\n";
			print $fh3 "$contig\t$value\t$contig\n"; #$plit[0]\n";
			#print $fh3 "$contig\t$value\t$split[0]\n";
			#print $fh3 "$contig\t$value\t$contig\n";
			push (@tags, $contig); #$split[0]);
		}
	}
}

@tags = uniq @tags;
foreach (@tags) {
	$_ =~ s/^\s+//;
	$_ =~ s/\s+$//;
	print $fh4 "$_\n";
}

close $fh4 || die "$!\n";
close $fh3 || die "$!\n";
