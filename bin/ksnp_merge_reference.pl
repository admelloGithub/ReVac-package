#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

=head1 NAME

ksnp_merge_reference.pl 

B<--input,-in,-i>
    single fasta file or list of fasta files

=head1  DESCRIPTION

Merge multiple reference genomes into one multifasta file for kSNP component.  Add "NN" if reference is made up of contigs/scaffolds

=head1  INPUT

fsa or list of fasta files

=head1  OUTPUT

A multi-fasta file that will be used as a reference genome in kSNP
    
=head1  CONTACT

    Kent Shefchek
    kshefchek@som.umaryland.edu

=cut

use strict;
use File::Temp;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

my $input;
my %options;
my $bool;
my $bool2 = 'false';

my $results = GetOptions(\%options,
			"input|in|i=s" => \$input
						 );

die "Must enter input file" unless defined $input;

open (my $fh, "<$input") || die "Cannot open input file";

if ($input =~ /.*\.list/) {
	my @list = <$fh>;
	foreach my $i (@list) {
		my $tmp = File::Temp->new(TEMPLATE => "tempXXXX",
			          DIR => "./",
			          SUFFIX => ".fsa"
			          );
		open (my $fh3, ">$tmp") || die "Cannot open temp file";		
		$i =~ s/\n//;
		open (my $fh2, "<$i") || die "Cannot open input in list file";
		my @fasta = <$fh2>;
		$bool = "true";
		foreach my $a (@fasta) {#Iterate and join contigs if needed
			if (($a =~ /^>/)&&($bool eq 'true')) {
				print $fh3 "$a";
				$bool = 'false';
			} elsif (($a =~ /^>/)&&($bool eq 'false')) {
				print $fh3 "NN";
				$bool2 = 'true';
			} else {
				$a =~ s/\n//;
				print $fh3 $a;
			}
		}
		print $fh3 "\n";
		close $fh3;
		#Now print with 60 bases per line

		open (my $fh4, "<$tmp") || die "Cannot open temp file";	
		
		my $base = fileparse($i, ".list");
		
		while (<$fh4>) {
			$_ =~ s/^>.*/>$base merged/ if ($bool2 eq 'true');
			$_ =~ s/(\w{60})/$&\n/g;
			$_ =~ s/^\s+$//g;
			print $_;
		}
	}
} elsif ($input =~ /.*\.f((a?st?)|n)?a/) { #Single fasta file
	my @list = <$fh>;
	my $tmp = File::Temp->new(TEMPLATE => "tempXXXX",
			          DIR => "./",
			          SUFFIX => ".fsa"
			          );
	open (my $fh3, ">$tmp") || die "Cannot open temp file";	
	$bool = "true";
	foreach my $a (@list) { #Iterate and join contigs if needed
		if (($a =~ /^>/)&&($bool eq 'true')) {
			print $fh3 "$a";
			$bool = 'false';
		} elsif (($a =~ /^>/)&&($bool eq 'false')) {
			print $fh3 "NN";
			$bool2 = 'true';
		} else {
			$a =~ s/\n//;
			print $fh3 $a;
		}
	}
	print $fh3 "\n";
	close $fh3;
	#Now print with 60 bases per line
	open (my $fh4, "<$tmp") || die "Cannot open temp file";	

	my $base = fileparse($input, qr/\.f((a?st?)|n)?a/);
	
	while (<$fh4>) {
		$_ =~ s/^>.*/>$base merged/ if ($bool2 eq 'true');
		$_ =~ s/^>.*/>$base/ if ($bool2 eq 'false');
		$_ =~ s/(\w{60})/$&\n/g;
		$_ =~ s/^\s+$//g;
		print $_;		
	}
} else { die "incorrect file format"; }










