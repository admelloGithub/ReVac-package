#!/usr/bin/perl -w
use strict;
use File::Basename;

my $python_path = $ARGV[0];
my $bcell_exec = $ARGV[1];
my $file = $ARGV[2];
my $base = basename($file);
$base =~ s/\.[a-z]+\.[a-z]+$//;

my @methods = qw(Emini Bepipred Parker Chou-Fasman Karplus-Schulz Kolaskar-Tongaonkar);

foreach (@methods) {
my $method = $_;
print "$base\t$python_path\t$bcell_exec\t$method\t$file\n";

}
