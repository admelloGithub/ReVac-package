#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
use strict;
use warnings;

#******************************************************************************
# *merge_uchime_ids.pl
# Author: james robert white, james.dna.white@gmail.com

# This function takes a list of files with chimera sequence reports from uchime
# and a list of files of clusters from mothur_unique_seqs as input.

# The output is a single concatenated list of ids that are all chimeric 
# sequences.
#******************************************************************************
use Getopt::Std;
use Data::Dumper;

use vars qw/$opt_n $opt_c $opt_o/;

getopts("n:c:o:");

my $usage = "Usage:  $0 \
                -n <list of name cluster files>\
                -c <list of files with representative chimera sequence ids>\
                -o <output file name>\
                \n";

die $usage unless defined $opt_n
              and defined $opt_c
              and defined $opt_o;

my $uniqueRepFileList = $opt_n;
my $chimeraFileList   = $opt_c;
my $outfile           = $opt_o;

my @uniqueclusterfiles = ();
my %clusters = ();

# load up list of unique read clusters
open (FLIST, $uniqueRepFileList) or die ("Could not open Unique read file list $uniqueRepFileList: $!");
while (my $file = <FLIST>) {
  chomp($file);
  push (@uniqueclusterfiles, $file);
}
close (FLIST);

# create a hash of all clusters
my $qiimeck = 0;
foreach my $file (@uniqueclusterfiles){
  open IN, $file or die "Could not open file $file!!\n";
  while(<IN>){
    chomp($_);
    my @A = split "\t", $_;
    if ($qiimeck == 0){
      if ($A[0] == 0){  # then we've got qiime formatted OTUs
        $qiimeck = 1;   
      }else{
        $qiimeck = -1;
      }
    }

    if ($qiimeck == -1){
      my @B = split ",", $A[1];
      for my $i (0 .. $#B){
        push @{$clusters{$B[0]}}, $B[$i];
      } # this cluster includes the representative in the hash value
    }else{
      for my $i (1 .. $#A){
        push @{$clusters{$A[0]}}, $A[$i];
      }
    }

  }
}

# now go through the list of chimera accnos and 
# print out all chimera ids using the clusters hash
my @chifiles = ();
open (FLIST, $chimeraFileList) or die ("Could not open Unique read file list $chimeraFileList: $!");
while (my $file = <FLIST>) {
  chomp($file);
  push (@chifiles, $file);
}
close (FLIST);


open OUT, ">$outfile" or die;
foreach my $file (@chifiles){ # for each chimera report file
  open IN, $file or die "Could not open file $file!!\n";
  while(<IN>){
    chomp($_);
    my @A = split "\t", $_;
    next if ($A[$#A] eq "N"); # b/c it's not detected as a chimera

    if ($qiimeck == 1){
      my @B = split /\//, $A[1];
      $A[1] = $B[0];
    }

    if (!defined($clusters{$A[1]})){
      die "Cannot locate representative sequence: $A[1]!!\n"; 
    }else{ #print that cluster out b/c theyre all chimeras
      print OUT "$A[1]\n";
      foreach my $v (@{$clusters{$A[1]}}){
        print OUT "$v\n"; 
      }
    } 
  }
  close IN;
}
close OUT;


