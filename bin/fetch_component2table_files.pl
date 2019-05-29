#!/usr/local/bin/perl -w 
use strict;
use Getopt::Std;
use File::Basename;
use List::MoreUtils qw(uniq);

#Example Usage:
#./fetch_component2table_files.pl -o /usr/local/projects/PNTHI/output_repository -p 15766 -l lipop,antigenic,mhc_class_i_all,
#																							immunogenicity,mhc_class_ii_all,tmhmm,autoimmunity,hmmpfam3,spaan,psortb,signalp,bcell_pred_all

getopts("o:p:l:");
our ($opt_o,$opt_p,$opt_l);

my @output_repo = split(/,/,$opt_o);
my $pipeline_id = $opt_p;
my @components = split(/,/,$opt_l);
my (@dirs,@paths);

foreach my $output_repo (@output_repo) {
	foreach (@components) {
		my $dp = `find $output_repo/$_/ -name "$pipeline_id\_*"`;
		if ($dp) {
		push (@dirs,$dp);
		} else {
			#my $dp = `find /usr/local/scratch/admello/output_repository/$_/ -name "$pipeline_id\_*"` || 
			warn "Directory may not have been found for $_\n";
			#push (@dirs,$dp);
		}
	}
}

my $dirs = join("\n",@dirs);
@dirs = uniq (split(/\n/,$dirs));
my (@hmms,%hmm,$odir);
foreach my $dir (@dirs) {
	if ($dir =~ /hmmpfam3/ && $dir !~ /\.old/ && $dir !~ /attributor/) {
		open (my $fh, "<", "$dir/hmmpfam32table.out") || die "$!:$dir/hmmpfam32table.out\n";
			my @lines = <$fh>;
		close $fh || die "$!\n";
		foreach (@lines) {
			my @split = split(/\t/,$_);
			$split[1] =~ s/\n//;
			$hmm{$split[0]} .= "$split[1];";
		}
	$odir = dirname($dir);
	$dir = "";	
	}
}

open (my $fh2, ">", "$odir/hmm_merge2table.out") || die "$!\n";
	#print $fh2 "Locus_tag\tHMM\n";
foreach (sort keys %hmm) {
	print $fh2 "$_\t$hmm{$_}\n";
}
close $fh2 || die "$!\n";

push (@dirs,"$odir/hmm_merge2table.out");

foreach (@dirs) {
	if ($_ =~ /\.old/) {next;}
	$_ =~ s/\s+$//;
	$_ =~ s/^\s+//;
	if ($_ ne "") {
		my $fp = `find $_ -maxdepth 1 -name *2table.out*`;
		if ($fp) {
		print "$fp";
		} else {
			warn "File not found for $_\n";
		}
	}
}

