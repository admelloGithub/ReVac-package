#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;

my $file = shift @ARGV;
my $outputfile = shift @ARGV;
open FILE,"$file" or die "Can't open $file";
my @files = <FILE>;
chomp @files;

#Multi-key sort order is not returning as expected with either LC_ALL=C or LC_COLLATE=C
$ENV{'LC_ALL'} = undef;
$ENV{'LC_COLLATE'} = undef;
$ENV{'LANG'} = 'en_US.UTF-8';
if(@files) {
    my $cmd = "sort -m ".join(' ',@ARGV)." ".join(' ',@files). " > $outputfile";
    print STDERR "Running $cmd\n";
    print `$cmd`;
}
