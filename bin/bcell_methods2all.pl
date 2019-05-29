#!/usr/bin/perl -w
use strict;
use File::Basename;
use List::MoreUtils qw(uniq);

my (%file,%table,%final,%peptides,@methods,%filebases);
my $methods = "Emini Bepipred Parker Chou-Fasman Karplus-Schulz Kolaskar-Tongaonkar Classical";
my $bcell_all_out_list = $ARGV[0];
my $outdir = $ARGV[1];

open (my $fh3, "<", $bcell_all_out_list) || die "$!\n";
my @lines = <$fh3>;
close $fh3 || die "$!\n";

foreach my $line (@lines) {
	my $base = basename($line);
	my @split = split (/\./,$base);
	$base = $split[0];
	$line =~ s/\s+$//;
	$filebases{$base} .= "$line,";
}


foreach my $base (keys %filebases) {

	my @files = split(/,/,$filebases{$base});
	my @seq_ids;
	#if (-e "$outdir/$base.bcell_pred_all") {next;}
	
		foreach (@files) {

			my $file = $_;
			open (my $fh, "<", $file) || die "$!:$file\n";
			my ($method,$seq_id);

			while (<$fh>) {
		
				if ($_ =~ /^>/) {$seq_id = $_; push(@seq_ids, $seq_id); next;}

				if ($_ eq "\n") {next;}

				if ($_ =~ /Position/) {
					my @split = split(/\s+/, $_);
					foreach (@split) {
						if ($methods =~ /$_/) {
							$method = $_;
							push(@methods, $_);
						}
					}
					next;
				}

				my $score;
				my @split = split(/\s+/, $_);
					foreach (@split) {
						if ($method =~ /Bepipred/) {
							$score = $split[$#split-1];
						} else {
							$score = $split[$#split];
						}
					}
				$table{$seq_id}{$split[0]}{$method} = "$score\t"; 
				
				if ($method eq "Parker") {
					$peptides{$seq_id}{$split[0]} = $split[$#split-1];
				}

			}

			close $fh || die "$!\n";

		}

	open (my $fh2, ">", "$outdir/$base.bcell_pred_all") || die ":$!\n";
		
		my @uniq = uniq @methods;
		my @uniq_seq_ids = uniq @seq_ids;
		foreach my $id (@uniq_seq_ids) {

			print $fh2 "$id"."Position ".join("\t",@uniq)."\n";

			my @positions = sort {$a <=> $b} keys %{$table{$id}}; 
			
			foreach my $pos (@positions) {

				if ($pos > 3 && $pos <$#positions-3) {

					if (!($table{$id}{$pos}{"Bepipred"})||!($table{$id}{$pos}{"Classical"})||!($table{$id}{$pos}{"Parker"})||!($table{$id}{$pos}{"Karplus-Schulz"})||!($table{$id}{$pos}{"Kolaskar-Tongaonkar"})||!($table{$id}{$pos}{"Emini"}))
					 { warn "$base\t$id\t$pos\t$table{$id}{$pos}{\"Bepipred\"}.$table{$id}{$pos}{\"Classical\"}.$table{$id}{$pos}{\"Parker\"}.$table{$id}{$pos}{\"Karplus-Schulz\"}.$table{$id}{$pos}{\"Kolaskar-Tongaonkar\"}.$table{$id}{$pos}{\"Emini\"}.\n"; }

					if ($table{$id}{$pos}{"Bepipred"}>=0.35 && 
						$table{$id}{$pos}{"Classical"}>=1.001 && 
						$table{$id}{$pos}{"Parker"}>=1.673 && 
						$table{$id}{$pos}{"Karplus-Schulz"}>=1 && 
						$table{$id}{$pos}{"Kolaskar-Tongaonkar"}>=1 &&
						$table{$id}{$pos}{"Emini"}>=1) {

						#print "$pos\t";
						my $start = $pos-3;
						my $end = $pos+3;
						print $fh2 "$start $end\t";	

						foreach my $key (@uniq) {
							
								if (!$table{$id}{$pos}{$key}) {
									print $fh2 "X.XXX\t\t";
								} else {
									print $fh2 "$table{$id}{$pos}{$key}\t"; 
								}
										
						}

						if ($peptides{$id}{$pos}) {
						print $fh2 "$peptides{$id}{$pos}\n";
						} else {
						print $fh2 "\n";
						}

					} 

				}

			}

		}	

	close $fh2 || die "$!\n";

}





