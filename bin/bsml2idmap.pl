#!/usr/bin/perl -w
use strict;

my $list = $ARGV[0];
open (my $fh1, "<", $list) || die "$!\n";
my @files = <$fh1>;
close $fh1;

foreach my $file (@files) {

	$file =~ s/\s+$//;
	`gzip -d $file`;
	$file =~ s/\.gz//;
	open (my $fh, "<", $file) || die "$!\n";

	my @lines = <$fh>;

	close $fh;
	`gzip $file`;

	#foreach my $i (0 .. $#lines) {
	#	if ($lines[$i] !~ /Feature-group/) {
	#		delete $lines[$i];
	#	}
	#}
	#
	#print scalar(@lines);
	#print $lines[1];

	foreach my $i (0 .. $#lines) {

		if ($lines[$i] eq "") {next;}
		if ($lines[$i] =~ /group-set/ && $lines[$i+5] =~ /polypeptide/) {
			my $lt = $lines[$i];
			my $pp = $lines[$i+5];
			$lt =~ s/^.*group-set="locus_tag_//;
			$lt =~ s/"\>\n$//;
			$pp =~ s/^.*featref="//;
			$pp =~ s/"\>\<\/Feature-group-member\>\n$//;
			print "$lt\t$pp\n";
			#print "$lines[$i]\n$lines[$i+5]";
		}
	}

}
