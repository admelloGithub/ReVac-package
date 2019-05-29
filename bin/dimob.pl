#!/usr/bin/perl
#Predict genomic islands using IslandPath-Dimob
use strict;
use warnings;


# changed for IGS environment
#use lib './dinuc_dimob';
use lib '/usr/local/packages/islandpath/dinuc_dimob';

use mobgene;
use genomicislands;
use Getopt::Long;

my $hmm_evalue= 100;

# changed for IGS environment
#my $hmmdb3 = './hmmpfam/Pfam_ls_mobgene_selected_June172005_HMMER3';
#my $hmmdb2 = './hmmpfam/Pfam_ls_mobgene_selected_June172005_HMMER2';
my $hmmdb3 = '/usr/local/packages/islandpath/hmmpfam/Pfam_ls_mobgene_selected_June172005_HMMER3';
my $hmmdb2 = '/usr/local/packages/islandpath/hmmpfam/Pfam_ls_mobgene_selected_June172005_HMMER2';

my $hmmer2=0;
GetOptions ("hmmer2" => \$hmmer2);

my $usage = "Usage:\n./dimob.pl [--hmmer2] <faa file> <ffn file> <ptt file>\nExample:\n./dimob.pl example/NC_000913.faa example/NC_000913.ffn example/NC_000913.ptt > islands.txt\n";

my ($faainput,$ffninput,$pttinput) = @ARGV;

unless(defined($faainput) && defined($ffninput) && defined($pttinput)){
    print $usage;
    exit;
}

my $hmminput = $faainput . ".mob";

#unless ( -e $hmminput ) { 											#############This is set up outside this script
    #Scan for mobility genes using hmmer (unless the file already exists)
#    if($hmmer2){
#        system("hmmpfam $hmmdb2 $faainput > $hmminput");
#    }else{
#        system("hmmscan $hmmdb3 $faainput > $hmminput");
#    }
#}

my $mob_list;

my $mobgenes = parse_hmmer( $hmminput, $hmm_evalue );

foreach(keys %$mobgenes){
   $mob_list->{$_}=1;
}

#get a list of mobility genes from ptt file based on keyword match
my $mobgene_ptt = parse_ptt($pttinput);

foreach(keys %$mobgene_ptt){
   $mob_list->{$_}=1;
}

#calculate the dinuc bias for each gene cluster of 6 genes
#input is a fasta file of ORF nucleotide sequences
my $dinuc_results = cal_dinuc($ffninput);
my @dinuc_values;
foreach my $val (@$dinuc_results) {
    push @dinuc_values, $val->{'DINUC_bias'};
}

#calculate the mean and std deviation of the dinuc values
my $mean = cal_mean( \@dinuc_values );
my $sd   = cal_stddev( \@dinuc_values );

#generate a list of dinuc islands with ffn fasta file def line as the hash key
my $gi_orfs = dinuc_islands( $dinuc_results, $mean, $sd, 8);

#convert the def line to gi numbers (the data structure is maintained)
my $dinuc_islands = defline2gi( $gi_orfs, $pttinput );

#check the dinuc islands against the mobility gene list
#any dinuc islands containing >=1 mobility gene are classified as
#dimob islands
my $dimob_islands = dimob_islands( $dinuc_islands, $mob_list );


foreach (@$dimob_islands) {

    #get the pids from the  for just the start and end genes
    my $start = $_->[0]{start};
    my $end = $_->[-1]{end};
    
    #print "$start\t$end\n";

	open (my $FH, "<", $pttinput) || die "Couldn't open ptt file for reading : $!\n";

		while (<$FH>) {
			if ($. > 3) {
				my @split = split (/\s+/,$_);
				my @coords = split (/\.\./,$split[0]);
				#$coords[0] =~ s/:[c]{0,}//;
				#print "$split[0]\t";#@coords\n"; 
					if ($coords[1] <= $end && $coords[0] > $start) {
						my $seqid = $split[0]; #join(" ",$split[1..4]);
						#$seqid =~ s/^>//;
						#$seqid =~ s/\.CDS\.[0-9]+\.1\|:[c]{0,1}[0-9]+-[0-9]+\n//;
						my $annot = join("_", @split[8..$#split]);
						$seqid =~ s/\.\./-/;
						my $id = `cat $faainput |grep "$seqid"`;
							if (!($id)) {
								my @split = split(/-/,$seqid);
								my $rev_coords = "$split[1]-$split[0]";
								$id = `cat $faainput |grep "$rev_coords"`;
							}
						my $id_coords = $id;
						$id_coords =~ s/.*\|//;
						$id_coords =~ s/\n//;
						#$id =~ s/\..*$//;
						$id =~ s/>//;
						$id =~ s/\n//;
						if (!($id)) { next; }
						print "$id\t$annot|$id_coords|GI:$start-$end\n";
					}
			}
		}

	close $FH || die "$!\n";	

}


#######################










