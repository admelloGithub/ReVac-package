#!/usr/bin/perl

## Function: Make a surface table with various features in different files
## Needs one tab-delimited file per feature, 2 columns: ORF \t Feature
## The -n option allows to include larger files for which we want to include
## the feature only if there at least another feature in the XXX columns after ORF

use Getopt::Std;
use File::Basename;

if ($#ARGV < 0) {
    print "\nNeeds at least 1 file: [-n XXX] file1 [file2 file3 ...]\n"
	."Files must be tab-delimited\n"
	."-n XXX indicates the lines to delete if XXX columns after ORF are all 'None'\n\n";
    exit;
}

getopts ('n:');

foreach $input (@ARGV) {
	$base = basename($input);
	$base =~ s/2table\.out.{0,}//;
	$table{head}{$input} = $base;
    open (INPUT, "$input") || die "Cannot open $input : $!\n";
    while(<INPUT>) {
	chomp;
	my @tabs = split(/\t/, $_);
		if (scalar(@tabs) > 2) {
			$orf = $tabs[0];
			$orf =~ s/\s+$//;
			shift @tabs;
			$feat = join("\t",@tabs);
		} else {
			($orf,$feat) = split(/\t/);
			$orf =~ s/\s+$//;
		}
		$table{$orf}{$input} = $feat;
    }
    close (INPUT) || die "Cannot close $input : $!\n";
}


if ($opt_n) {
    foreach $key1 (keys %table) {
	foreach $key2 (@ARGV) {
	    $ct ++;
	    $string = $table{$key1}{$key2};
	    unless ($string) {
		$nb ++;
		if ($ct == $opt_n and $nb == $opt_n) {
		    $junk{$key1} = "Delete";
		}
	    }
	}
	$ct = 0;
	$nb = 0;
    }
}

print "Locus_Tag\t";
foreach (sort @ARGV) {
	print "$table{head}{$_}\t";
}
print "\n";

foreach $key1 (keys %table) {
	if ($key1 eq "head") {next;}
    unless (defined $junk{$key1}) {
	print "$key1\t";
	foreach $key2 (sort @ARGV) {
	    $string = $table{$key1}{$key2};
	    if ($string) {
		print "$string\t";
	    }
	    else {
		print "None\t";
	    }
	}
	print "\n";
    }

}

__END__

