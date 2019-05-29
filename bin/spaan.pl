#!/usr/bin/perl

use strict;
use warnings;
use Cwd;
use File::Copy;
use File::Basename;

my $pwd = cwd();
my $input = $ARGV[0]; #Accepts cmd line split_multifasta list file as input
my $output = $ARGV[1]; #Path to output folder with some file name
#my @file;
#system("rm -f $output/spaan.output");
#system("> $output/spaan.output");

if (-e $ARGV[0]) {
my @file;

if ($ARGV[2]) {
	print "Please enter a single input file followed by output file destination and name.\n";
	exit(1);
}


if ($input =~ /.list$/) {

	open (my $fh, "<", "$input") 
		or die "Failed to open bsml2fasta.multi list:$!\n";

	@file = <$fh>;

	close $fh
 	or die "Failed to close bsml2fasta.multi list:$!\n";
	
	foreach(@file) {
	$_ =~ s/^\s+//;
	$_ =~ s/\s+$//;
	}
	
}

foreach(@file) {

    $input = $_;
    my $basename = basename($_);

if ($input =~ /.fsa$/) {
		copy("$input","/usr/local/projects/PNTHI/tools/SPAAN/query.dat") or die "Copy failed: $!";  #SPAAN runs from local input query.dat
}

else {

	print "Please enter a valid fasta formatted .fsa file\n";
	exit(1);
}


if (copy("$input","/usr/local/projects/PNTHI/tools/SPAAN/query.dat") == 1) {

	chdir("/usr/local/projects/PNTHI/tools/SPAAN") or die "$!";

        system("/usr/local/projects/PNTHI/tools/SPAAN/askquery");

	copy("/usr/local/projects/PNTHI/tools/SPAAN/query.out","$output/$basename.spaan") or die "Copy failed: $!";

	system("rm query.out");
	system("rm query.dat");

	#print "Check /home/admello/Desktop/spaan/spaan.output \n";

	chdir("$pwd") or die "$!";

}

else {
	print "Invalid file!\n";
	exit(1);
}

}
}
################################
