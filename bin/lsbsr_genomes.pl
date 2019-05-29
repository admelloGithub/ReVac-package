#!/usr/bin/perl -w
use strict;
use File::Basename;
use List::MoreUtils qw(uniq);
##Appends .fasta extension to files and copies them into a given directory.
##Joins multiple contig fastas into single strain fastas.

my $list = $ARGV[0];
my $outdir =$ARGV[1];

if (!(-e $outdir)) {
	system ("mkdir -p -m 777 $outdir");
} 

open (my $fh, "<", $list) || die "$!\n";

my @files = <$fh>;

close $fh || die "$!\n";

my @keys;

foreach (@files) {
    $_ =~ s/\s+$//;
    if (-e $_) {
        my $file = $_;
        $file = basename($file);
        $file =~ s/\..*$//;
		my @k;
		#if ($file =~ /0{3,}/ && $file =~ /_/) {  @k  = split(/0{3,}/,$file); }
		#els
		if ($file =~ /0{45,}/) {	@k  = split(/0{3,}/,$file); } else { $k[0] = $file; }
		#elsif ($file =~ /_/) {	@k = split(/_/,$file); }
		#if ($file !~ /0{3,}|_/) { push (@keys, $file); }
		if ($file !~ /0{3,}|_/) {die "File names aren't in the required format $file\n";	@k = split(/\./,$file); }
		push (@keys, $k[0]);
	}
}

my @newfiles = uniq (@keys);
foreach (@newfiles) {
print $_."\n"; #scalar(@newfiles);
}
foreach my $nf (@newfiles) {	
	foreach (uniq @files) {
		$_ =~ s/\s+$//;
		if (-e $_) {
			my $file = $_;
			$file = basename($file);
			#$file =~ s/\..*$//;
			if ($file =~ /$nf/) {
			#$file .= ".fasta";
			system ("cat $_ >>$outdir/$nf.fasta");
			#system ("cp $_ $outdir/$file");
			}
		} else {
			die "List file may not contain actual files.\n";
		}
	}
}

exit;
