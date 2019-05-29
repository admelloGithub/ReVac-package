#!/usr/bin/perl -w
use strict;
use File::Basename;

my $gbk_file = $ARGV[0];
my $locus_tag_prefix = $ARGV[1];
my (@seqs,%genes,%annot,%loci,%coords,%contig,$contig);
my $ct = 0;
my $j = 0;
my ($title,$strain,$isolate,$date);

open (my $fh, "<", $gbk_file) || die "Couldn't open genbank file:$!\n";

while (<$fh>) {

	if ($_ =~ /LOCUS/) {
		my @split = split (/\s+/,$_);
		$contig = $split[1];
	}

	
		if ($_ =~ /\/organism=/) {
		$title = $_;
		$title =~ s/.*\/organism=//;
		$title =~ s/ /_/g;
		$title =~ s/"//g;
		$title =~ s/\s+$//;
		}
		if ($_ =~ /\/strain=/) {
		$_ =~ s/.*\/strain=//;
		$_ =~ s/ /_/g;
		$_ =~ s/"//g;
		$_ =~ s/\s+$//;
		#$strain ="";
		$strain = "_$_" if ($title !~ /$_/);
		}
		if ($_ =~ /\/isolate=/) {
		$_ =~ s/.*\/isolate=//;
		$_ =~ s/ /_/g;
		$_ =~ s/"//g;
		$_ =~ s/\s+$//;
		#$isolate = "";
		$isolate = "_$_" if ($title !~ /$_/);
		}
		if ($_ =~ /\/collection_date=/) {
		$_ =~ s/.*\/collection_date=//;
		$_ =~ s/"//g;
		$_ =~ s/^.*-//;
		$_ =~ s/\//or/g;
		$_ =~ s/\s+$//;
		#$date = "";
		$date = "_$_";
		}


	if ($_ =~ /^\s+CDS\s{2,}/ && $_ !~ /::/) {
		$contig{$contig} .= $_;		
			while (<$fh>) { 
 
				if ($_ !~ /ORIGIN/ && $_ !~ /CONTIG/ && $_ !~ /[ ]{2,}t[m]{0,}RNA[ ]{2,}/ && $_ !~ /^\s+gene[ ]{2,}/) {    ##Ensures only stores CDS sections of gbk
					$contig{$contig} .= "$_";
				} else { last; }
				if ($_ =~ /translation=/) { $ct = 1; }		#Ensures translation is parsed in event previous feature ends as like a protein terminal end
				if ($_ =~ /^\s+[ARNDBCEQZGHILKMFPSTWYVX]{0,}"$/ || $_ =~ /^\s+\/translation="[ARNDBCEQZGHILKMFPSTWYVX]{0,}"$/) {
					if ($ct ==1) {
						$ct = 0;
						last;
					}
				}

			}
	}

}

close $fh ||die "$!\n";

foreach my $contig (sort keys %contig) {

	my $lines = $contig{$contig};
	$lines =~ s/^\s+CDS/||/g;
	$lines =~ s/\n\s+CDS\s+/||/g;
	$lines =~ s/\n//g;
	$lines =~ s/[ ]{2,}//g;
	$lines =~ s/\t//g;
	$lines =~ s/"//g;

	my @lines = split(/\|\|/,$lines);
	shift @lines;

	my @aaa;
	foreach (@lines) {

		if ($_ =~ /\//) {

		my @feats = split(/\//,$_);
		my $seq = $feats[$#feats];
		if ($seq =~ /translation=/) {
			$seq =~ s/translation=//;
		} else {
			push(@aaa,$seq);
			warn "Encountered incomplete CDS capture: All sequences might not be parsed\n$seq\n";
			next;
		}

		$feats[0] =~ s/complement/c/ if ($feats[0] && $feats[0] =~ /complement/); 
		
		if ($coords{$seq}) { 
			$j++;
			$seq .= "$j";		##In the event that proteins have two locus tags. Prevents swapping values in keys.
		}
			$coords{$seq} = "$contig\_$feats[0]";

				foreach (@feats) {
					if ($_ =~ /^gene=/) {
						$_ =~ s/gene=//;
						$genes{$seq} = $_;
					}
					if ($_ =~ /^locus_tag=/) {
						$_ =~ s/locus_tag=//;		
						$loci{$seq} = "N$_";
					}
					if ($_ =~ /^product=/) {
						$_ =~ s/product=//;
						$_ =~ s/\s+/_/g;
						$annot{$seq} = $_;
					}
				}
			
		push (@seqs, $seq);
		
		} else { next; }#die "Processing error= Delimiter \"\\\" not found in processed gbk line:$_\n";}

	}

}
#print scalar(@seqs)."\n";
#my $name = basename($ARGV[0]);
#$name =~ s/\.gff.*$//;
#print $name."\n";
#foreach (@aaa) {
#print ">$_\n";
#}

foreach (@seqs) {
	my $s_id = $_; 						###This is the entire protein.
	$s_id =~ s/[0-9]{0,}$//;			###Removes added numbers from ends of protein sequence
	print "$loci{$_}\t";
	print "$title";
    print "$strain" if ($strain);
    print "$isolate" if ($isolate);
    print "$date" if ($date);
	print "\t$coords{$_}\t";
	print length($s_id)."\t";
	#print calc_mass($s_id)."\t";
	print mw_calc($s_id)."/".iep_calc($s_id)."\t";
	my $seq = sprintf "%-60s", $annot{$_};
	print "$seq\t";
	if ($genes{$_}) {print "$genes{$_}\t";} else {print "None\t";}
	#print "$name";
	print "\n";
}

############Calculates protein molecular weight##########
####BELOW SUBROUTINES ACQUIRED FROM STACK OVERFLOW####
sub calc_mass {
    my $a = shift;
    my @a = ();
    my $x = length $a;
    @a = split q{}, $a;
    my $b = 0;
    my %data = (
        A=>71.09,  R=>16.19,  D=>114.11,  N=>115.09,
        C=>103.15,  E=>129.12,  Q=>128.14,  G=>57.05,
        H=>137.14,  I=>113.16,  L=>113.16,  K=>128.17,
        M=>131.19,  F=>147.18,  P=>97.12,  S=>87.08,
        T=>101.11,  W=>186.12,  Y=>163.18,  V=>99.14, X=>0.00
    );
    for my $i( @a ) {
        $b += $data{$i};
    }
    my $c = ($b - (18 * ($x - 1)))/100;
    my $KDa = sprintf "%.2f", $c;
	$KDa = $KDa."KDa";
	$KDa = sprintf "%-10s", $KDa;
	return $KDa;
}

  sub mw_calc {
     my $peptide    = shift;
     my $kd       = shift || 0;;
      my $mass=0;
     my $index;
     my %aas;
 #if ($kd){ #Kyte-Doolittle
     %aas = (
 A   => 71.0788,
 B => 114.6686,   # Asx   Aspartic acid or Asparagine 
 C   => 103.1388,
 D   => 115.0886,
 E   => 129.1155,
 F   => 147.1766,
 G   => 57.0519,
 H   => 137.1411,
 I   => 113.1594,
 K   => 128.1741,
 L   => 113.1594,
 M   => 131.1926,
 N   => 114.1038,
 O   => 237.3018,
 P   => 97.1167,
 Q   => 128.1307,
 R   => 156.1875,
 S   => 87.0782,
 T   => 101.1051,
 U   => 150.0388,
 V   => 99.1326,
 W   => 186.2132,
 X => 111.1138, # Xaa   Any amino acid
 Y   => 163.176,
 Z => 128.7531    #Glx   Glutamine or Glutamic acid
          );
#  }
  my $len = length($peptide);
 
     for ($index = 0; $index < $len; $index++) {
 
     my $letter = substr($peptide, $index, 1);
     $letter=uc($letter);
 
     $mass += $aas{$letter};
 
     }
 $mass=$mass+18.015;
 $mass+=0.0005;
 $mass*=1000;
 $mass=~s/\..+//;
 $mass=$mass/1000000;
 $mass = substr($mass, 0, 5);
     return $mass;
}

sub iep_calc {

# Total number N of the amino acids Lys(K),Arg(R),His(H),Asp(D),Glu(E),Cys(C)and Tyr(Y) within the sequence.

my $N_K=0;
my $N_R=0;
my $N_H=0;
my $N_D=0;
my $N_E=0;
my $N_C=0;
my $N_Y=0;


my $sequence=shift;
while ($sequence =~ /K/ig){$N_K++}
while ($sequence =~ /R/ig){$N_R++}
while ($sequence =~ /H/ig){$N_H++}
while ($sequence =~ /D/ig){$N_D++}
while ($sequence =~ /E/ig){$N_E++}
while ($sequence =~ /C/ig){$N_C++}
while ($sequence =~ /Y/ig){$N_Y++}

# print fileout "K = $N_K\n";
# print fileout "R = $N_R\n";
# print fileout "H = $N_H\n";
# print fileout "D = $N_D\n";
# print fileout "E = $N_E\n";
# print fileout "C = $N_C\n";
# print fileout "Y = $N_Y\n\n";

# K values.

my $K_lys = 10**(-10.00);
my $K_arg = 10**(-12.00);
my $K_his = 10**(-5.98);
my $K_asp = 10**(-4.05);
my $K_glu = 10**(-4.45);
my $K_cys = 10**(-9.00);
my $K_tyr = 10**(-10.00);


my %pK_NT = (

A => 7.59,
C => 7.50,
D => 7.50,
E => 7.70,
F => 7.50,
G => 7.50,
H => 7.50,
I => 7.50,
K => 7.50,
L => 7.50,
M => 7.00,
N => 7.50,
P => 8.36,
Q => 7.50,
R => 7.50,
S => 6.93,
T => 6.82,
V => 7.44,
W => 7.50,
Y => 7.50
#B =>
#O =>
#U =>
#X =>
#Z =>
);

        my $NT_letter = substr ($sequence,0,1);
        $NT_letter=uc($NT_letter);
        my $K_NT = 0;
        $K_NT=10**(-$pK_NT{$NT_letter});

my %pK_CT = (

A => 3.55,
C => 3.55,
D => 4.55,
E => 4.75,
F => 3.55,
G => 3.55,
H => 3.55,
I => 3.55,
K => 3.55,
L => 3.55,
M => 3.55,
N => 3.55,
P => 3.55,
Q => 3.55,
R => 3.55,
S => 3.55,
T => 3.55,
V => 3.55,
W => 3.55,
Y => 3.55
);

        my $CT_letter = substr ($sequence,-1,1);
        $CT_letter=uc($CT_letter);
        my $K_CT = 10**(-$pK_CT{$CT_letter});
my $pH;
for ($pH=0;$pH<=14;$pH+=0.01){

        my $C_H = 10**(-$pH);

        my $CRp = $N_K*$C_H/($C_H+$K_lys) + $N_R*$C_H/($C_H+$K_arg) + $N_H*$C_H/($C_H+$K_his) + $C_H/($C_H+$K_NT);

        my $CRn = $N_D*$K_asp/($C_H+$K_asp) + $N_E*$K_glu/($C_H+$K_glu) + $N_C*$K_cys/($C_H+$K_cys) + $N_Y*$K_tyr/($C_H+$K_tyr) + $K_CT/($C_H+$K_CT);

        my $R = $CRp -$CRn;

        if (($pH==14) and ($R>=0))
        {
                return 14;
        }

        if ($R < 0){

        $pH -= 0.01;

$pH+=0.005;
$pH*=100;
$pH=~s/\..+//;
$pH/=100;

        return $pH;

                }
        }

}

