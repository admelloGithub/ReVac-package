#!/usr/bin/perl -w

# NERVE, copyright (C) 2006, author Sandro Vivona.

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, US

# author contact: sandro@bio.unipd.it


#useage ----   ./autoimmunity.pl -p1=records/Hemoglobin.fsa -o=/home/admello/Desktop/autoimmunity.output -db=/usr/local/projects/PNTHI/databases/



$orainizio = qx{date +%H:%M:%S};

print "
###################################################################################\n
                      statrting autoimmunity.pl at: $orainizio
###################################################################################\n
";



######################################################### ACQUISIZIONE PARAMETRI ############################################################

#settiamo i parametri di default: numero massimo di mismatch & substitution consentiti x ogni finestra di lunghezza voluta (lunghezza minima dei peptidi)
$min_length = 9;
$substitution = 10;
$mismatch = 10;



use File::Basename;




chomp($root=qx{pwd});
$Query = $1 if ($root =~ /\/([^\/]+)$/);
$Query =~ s/[ -\."']/_/g;

#non dovrebbe servire mai...pwd dovrebbe aver fatto il suo dovere!giusto x sicurezza!
if (!$root){
	print "\nAttention, I need your help to know your work directory address\n";
	print "For example /home/sandro/mypathogen\n";
	print "Please, type it now:";
	while (!$root){
		chomp($root = <STDIN>);
		if (!(-d "$root")){
			$root = "";
			print "\nSorry, this directory doesn't exist\n";
			print "Try again:";
		}
		else{
			print "\nOk, the address is correct!\n\n";
		}
	}
}

#passiamo al programma (come argomento) gli altri parametri dell'utente
foreach $arg (@ARGV){
	if ($arg =~ /^-proteome1=/ || $arg =~ /^-p1=/){
		@proteome1 = split ('=',$arg);
		$proteome1 = $proteome1[1]; #if ($proteome1[1] =~ /^\//);
		#$proteome1 = "$root/".$proteome1[1] if ($proteome1[1] !~ /^\//);
	}
#	elsif ($arg =~ /^-mysql=/ || $arg =~ /^-m=/){							###Modified for ergatis environment.
#		@mysql = split ('=',$arg);
#		($host,$database,$user,$password) = split (',',$mysql[1]);
#	}
	elsif ($arg =~ /^-minlength=/ || $arg =~ /^-ml=/){
		@min_length = split ('=',$arg);
		$min_length = $min_length[1];
	}
	elsif ($arg =~ /^-substitution=/ || $arg =~ /^-ss=/){
		@substitution = split ('=',$arg);
		$substitution = $substitution[1];
	}
	elsif ($arg =~ /^-mismatch=/ || $arg =~ /^-mm=/){
		@mismatch = split ('=',$arg);
		$mismatch = $mismatch[1];
	}
	elsif ($arg =~ /^-output_dir=/ || $arg =~ /^-o=/){
		@output_dir = split ('=',$arg);
		$output_dir = $output_dir[1];
	}
	elsif ($arg =~ /^-database=/ || $arg =~ /^-db=/){
                 @humandb_dir = split ('=',$arg);
                 $humandb_dir = $humandb_dir[1];
        }
	elsif ($arg =~ /^-blast_path=/ || $arg =~ /^-bp=/){
                 @blast_path = split ('=',$arg);
                 $blast_path = $blast_path[1];
        }

}


print "The following parameters have been set for autoimmunity.pl\n";
print "	root		=	$root\n";
print "	proteome1	=	$proteome1\n";
print "	min_len		=	$min_length\n";
print "	substitution	=	$substitution\n";
print "	param_mm	=	$mismatch\n\n";

$filename=basename($proteome1);

##################################


print "\n\n...extracting human entries from mhcpep...writing hladrpep...\n\n";
open (MHCPEP,"$humandb_dir/mhcpep.txt");
open (HLADRPEP,"+>$humandb_dir/hladrpep.txt");
while($rigaMHCPEP=<MHCPEP>){
        $reportMHCPEP.=$rigaMHCPEP;
}

@entries = split ('>',$reportMHCPEP);
foreach $e (@entries){
        print HLADRPEP ">$e";# if($e =~ /^HUM/);
}
if(-z "$humandb_dir/hladrpep.txt"){
        warn "\n$!\nParsing mhcpep human entries failed.\n";
       # print "Screenings versus HLA ligand database will not be possible.Hit any key to continue or type 'q' to quit the work\n";
        #chomp($risposta04 = <STDIN>);
        #die if $risposta04 eq "q";
}

###################################

#################################################### BLAST PATOGENO CONTRO H_sapiens ########################################################
unless (-e "$output_dir/blast") {
	mkdir("$output_dir/blast",0755)||warn "$!";
}
#opendir(INP,"$root/records")||die "cannot open la directory records:$!";
#while($ID = readdir INP){
#	$allID++ if($ID!~ /^\./);
#}
#open(INP,"$root/records")||die "cannot open la directory records:$!";
#while($ID=readdir INP){
#	if($ID!~ /^\./){
		
	#eval("formatdb -i $humandb_dir/H_sapiens.fasta");
		#print "blastall -p blastp -F N -i $proteome1 -d $humandb_dir/H_sapiens.fasta -o $output_dir/blast/$filename.blast \n";
		`$blast_path -p blastp -F N -i $proteome1 -d $humandb_dir/H_sapiens.fasta -o $output_dir/blast/$filename.blast`;
		
		$c++;
		print "Proteins in $filename just compared Vs Homo sapiens\n";#\t$c/$allID\n";
#	}
#}


############################################### CONNESSIONE AL DATABASE MYSQL E AGGIUNTA COLONNE ########################################



#use DBD::mysql;
#$db=Mysql->connect("$host","$database","$user","$password") || die "Fallita la connessione a Mysql!:$!\n";

#$db->query("alter table $Query add column length int after ID") || warn "Mysql column length creation failed:$!";

#$db->query("alter table $Query add column N_HomolPept2H_sapiens_l$min_length\_s$substitution\_m$mismatch smallint") || warn "Creation of Mysql column N_HomolPept2H_sapiens_l$min_length\_s$substitution\_m$mismatch failed:$!";

#$db->query("alter table $Query add column HomolPept2H_sapiens_l$min_length\_s$substitution\_m$mismatch text") || warn "Creation of Mysql column HomolPept2H_sapiens_l$min_length\_s$substitution\_m$mismatch failed:$!";

#$db->query("alter table $Query add column N_AA2_H_sapiens int") || warn "Mysql column N_AA2_H_sapiens creation failed:$!";


############################################# PREPARAZIONE DEL SET DI HLA-DR PER IL CONTROLLO ###########################################

open (H,"$humandb_dir/hladrpep.txt")||warn "could not open hladrpep.txt; screening to HLA-DR database will not be possible!";
while (<H>){
	push (@HLADR,$1) if(/SEQUENCE: (.+)\*\#/);
}

################################################## LETTURA ED ESTRAZIONE RISULTATI BLAST ################################################



#opendir (RIS,"$output_dir/blast/")||die "cannot open blast:$!";
#while($ID=readdir RIS){

#	next if($ID =~ /^\./);
	
	#$runningID++;
	#print "autoimmunity.pl processing $ID	$runningID/$allID\n";
	
	mkdir("$output_dir/HomolPept2H_sapiens_l$min_length\_s$substitution\_m$mismatch",0755);
	open (P,"+>$output_dir/HomolPept2H_sapiens_l$min_length\_s$substitution\_m$mismatch/$filename");

	print P "\n";
	print P "#******************************************************************************************#\n";
	print P "# +--------------------------------------------------------------------------------------+ #\n";
	print P "# | ************************************************************************************ | #\n";
	print P "# | #                                                                                  # | #\n";
	print P "# | #    This file, generated by the script autoimmunity.pl (author Sandro Vivona),    # | #\n";
	print P "# | #    reports the homologue peptides shared with Homo Sapiens. If you wish to       # | #\n";
	print P "# | #    design a subunit vaccine, you'd better discard them as potential cause of     # | #\n";
	print P "# | #    low immunogenicity or autoimmunogenicity problems. Theese peptides are        # | #\n";
	print P "# | #    consistent with the setting values you gave to autoimmunity.pl about number   # | #\n";
	print P "# | #    of mismatches & substitutions allowed for any amino acid window of desired    # | #\n";
	print P "# | #    length (the minimal length of peptides).                                      # | #\n";
	print P "# | #                                                                                  # | #\n";
	print P "# | #    Data are shown according to the underlying convetion                          # | #\n";
	print P "# | #                                                                                  # | #\n";
	print P "# | #            >human entry description                                              # | #\n";
	print P "# | #            ^pathogen peptide (Query)^                                            # | #\n";
	print P "# | #            *ncbiblast alignment output*                                          # | #\n";
	print P "# | #            #human peptide (Sbjct)#                                               # | #\n";
	print P "# | #                                                                                  # | #\n";
	print P "# | ************************************************************************************ | #\n";
	print P "# +--------------------------------------------------------------------------------------+ #\n";
	print P "#******************************************************************************************#\n";
	print P "\n";
	
	
	print P "\n";
	print P "############################################################################################\n";
	print P " See below the homologue peptides (not more than $mismatch mismatch and $substitution substitutions\n";
	print P " every $min_length aminoacid stretch) found comparing $filename proteins versus Homo Sapiens \n";
	print P "############################################################################################\n"; 
	
	
	#svuoto variabili che ci servono vergini per ogni ID
	$peptidi="\nno peptides responding to your parametres\n\n";
	$report= "";
	$c = 0;
	@peptidixAC = ();
 	@peptidixACnr = ();
	%freq_pep = ();
	$NpeptidixACnr = 0;
	@peptidixACnrTABfreq = ();
	$peptidixACnrTABfreq = "no peptides";
	$N_AA = 0;
	$ID_Length = 0;
	
	
	#comincio ad analizzare i file di blast
	open (R,"$output_dir/blast/$filename.blast")||die "$!";
	while(<R>){
		$report.=$_;
	}
	
	$report =~ /\((\d+) letters\)/;
	$ID_Length = $1;
#	$db->query("update $Query set length=$ID_Length where ID='$ID'")||warn "Insertion in mysql column length failed for protein $ID:$!";	


	if($report =~ /\*\*\*\*\* No hits found \*\*\*\*\*\*/){
		print P $peptidi;
		
#		$db->query("update $Query set N_HomolPept2H_sapiens_l$min_length\_s$substitution\_m$mismatch=$NpeptidixACnr where ID='$ID'")||warn "Insertion in mysql column N_HomolPept2H_sapiens_l$min_length\_s$substitution\_m$mismatch failed for protein $ID:$!";
		
#		$db->query("update $Query set HomolPept2H_sapiens_l$min_length\_s$substitution\_m$mismatch='$peptidixACnrTABfreq' where ID='$ID'")||warn "Insertion in mysql column HomolPept2H_sapiens_l$min_length\_s$substitution\_m$mismatch failed for protein $ID:$!";
		
#		$db->query("update $Query set N_AA2_H_sapiens=$N_AA where ID='$ID'")||warn "Insertion in mysql column N_AA2_H_sapiens failed for protein $ID:$!";
		
#		next;
	}
###########################################Tying IDs to blasted Human_genes

	
my @lines = split (/Query=/,$report); 
shift (@lines); 
foreach (@lines) {			##Closes at end of below foreach $align loop
my @split = split(/\n/, $_);
my $id = $split[0];			
	#######################################################################		
	@allineamenti = split(/^>/m,$_); #$report);
	shift (@allineamenti);
	
	foreach $align (@allineamenti){
	
		#riconoscimento AC&descriz della entry estratta
# 		$align =~ /\(([A-Z0-9]{6})\)/;
# 		$H_sapiensAC="$1";
		$align =~ /^(.+)\n/;
		$H_sapiensDescr="$1";

		#svuoto variabili che servono vergini ad ogni $align
		%H_sapienspeptidisorg = ();
		%Querypeptidisorg = ();
		
		#riconoscimento,pulitura e concatenamento in un unica stringona x le righe di allineamento condivise (così pure x le righe di Query e Sbjct)
		@parti_a = split (' Score =',$align);
		shift(@parti_a);
		foreach $parte (@parti_a){

			#svuoto variabili che servono vergini ad ogni $parte
			$stringona1 = "";
			$stringona = "";
			$stringona2 = "";
			$parte =~ /\n(Query: (\d+) +)\w/;
			$Nspaziiniziali = length($1);
			$Query_start_pos = $2;
			
			$parte =~ /\nSbjct: (\d+) +\w/;
			$H_sapiens_start_pos = $1;
			
			@righe = split('\n',$parte);
			foreach $r (@righe){
				if ($r =~ /(^ {$Nspaziiniziali}) *[A-Z\+ ]*/ && $r !~ /[a-z)(1-9,.:;]/){
					$r =~ s/$1//;
					$stringona.=$r;
				}
				elsif ($r =~ /^Query: \d+ +([\-A-Z]+) +\d+/){
					$stringona1.=$1;
				}
				elsif ($r =~ /^Sbjct: \d+ +([\-A-Z]+) +\d+/){
					$stringona2.=$1;
				}
			}
		
			#verifichiamo che i concatenamenti siano avvenuti correttamente
			$lstringona1 = length($stringona1);
			$lstringona = length($stringona);
			$lstringona2 = length($stringona2);
			print "attention an error occurred: $lstringona1 != $lstringona != $lstringona2. Extraction of source peptides from $H_sapiensDescr may not be precise\n" if ($lstringona1 != $lstringona || $lstringona != $lstringona2);
		
			#analisi della stringona di allineamento mediante finestra di scorrimento a lunghezza fissa (cuore del programma!!) ed estrazione dei peptidi e dei corrispondenti "peptidisorgenti"
			$z = 0;
			$stato = 0;
			for($h = 0; $h <= $lstringona-$min_length; $h++){
				
				$s = substr($stringona,$h,$min_length);
				@mm = ($s =~ /( )/g);
				@plus = ($s =~ /(\+)/g);
				
				if(@mm > $mismatch || @plus > $substitution){
					$i = $h-$z+$min_length-1;
					if($stato == 1){
						$peptide = substr ($stringona,$z,$i);
# 						$N_AA = $N_AA + length($peptide);
						$peptide = "\t*$peptide*";
						
						$Querypeptidesorg = substr ($stringona1,$z,$i);
						&contagap($stringona1,$z);$Query_pos = $Query_start_pos+$z-$gap;
						$Querypeptidesorg = "$Query_pos\t^$Querypeptidesorg^";
						$Querypeptidisorg{$peptide} = $Querypeptidesorg;push (@peptidixAC,$Querypeptidesorg);
						
						$H_sapienspeptidesorg = substr ($stringona2,$z,$i);
						&contagap($stringona2,$z);$H_sapiens_pos = $H_sapiens_start_pos+$z-$gap;
						$H_sapienspeptidesorg = "$H_sapiens_pos\t#$H_sapienspeptidesorg#";
						$H_sapienspeptidisorg{$peptide} = $H_sapienspeptidesorg;
					}
					$stato = 0;
				}
				else{
					$z = $h if($stato == 0);
					$stato = 1;
					if($h == $lstringona-$min_length){
						$i = $h-$z+$min_length-1 if (!$i);
						
						$peptide = substr ($stringona,$z,$i+1);
# 						$N_AA = $N_AA + length($peptide);
						$peptide = "\t*$peptide*";
						
						$Querypeptidesorg = substr ($stringona1,$z,$i+1);
						&contagap($stringona1,$z);$Query_pos = $Query_start_pos+$z-$gap;
						$Querypeptidesorg = "$Query_pos\t^$Querypeptidesorg^";
						$Querypeptidisorg{$peptide} = $Querypeptidesorg; push (@peptidixAC,$Querypeptidesorg);
						
						$H_sapienspeptidesorg = substr ($stringona2,$z,$i+1);
						&contagap($stringona2,$z);$H_sapiens_pos = $H_sapiens_start_pos+$z-$gap;
						$H_sapienspeptidesorg = "$H_sapiens_pos\t#$H_sapienspeptidesorg#";
						$H_sapienspeptidisorg{$peptide} = $H_sapienspeptidesorg;
					}
				}
			}
		}


		#scrivo nel file i peptidi selezionati e creo arrays associativi x stabilire frequenze da scrivere alla fine in fondo al file & x eliminar ridondanze (nr=non redundant) in vista dell'inserimento in tabella mysql


		if (%Querypeptidisorg){
			
			print P "\n>$H_sapiensDescr\n\n";
			#my $key = ">$H_sapiensDescr";
			#$key =~ s/\s+$//;
			print P "$id\n";  #"$seqs{$key}\n";
			foreach $p (sort { $Querypeptidisorg{$a} cmp $Querypeptidisorg{$b} }keys(%Querypeptidisorg)){
				print P "$Querypeptidisorg{$p}\n$p\n$H_sapienspeptidisorg{$p}\n\n";
			}

			$freq_H_sapiens{$H_sapiensDescr}++;
		
			foreach $pep (@peptidixAC){
				$freq_pep{$pep}++;
			}
			
			$c++;
			
		}
	}
	
	print P $peptidi if ($c == 0);

}
############################ PREPARAZIONE PEPTIDI X MYSQL (COMPRESA SCANSIONE CONTRO HLAD-DR DA MHCPEP) #################################

	@peptidixACnr = keys(%freq_pep);
	if (@peptidixACnr){
		$NpeptidixACnr = @peptidixACnr;
		foreach $pnr (@peptidixACnr){
			
			$contro_hladrpep = "negative";
			
			#tolgo momentaneamente cappelletti delimitatori e numeri-posizione x fare la scansione, ma prima conservo $pnr
			$backup = $pnr;
			$pnr =~ s/^\d+\t\^// || print "fallita rimozione cappelletto iniziale da $pnr\n";
			$pnr =~ s/\^$// || print "fallita rimozione cappelletto finale da $pnr\n";
			
			foreach $H (@HLADR){
				print "found a peptide positive to HLA_DR db\n" if ($contro_hladrpep eq "positive");
				last if ($contro_hladrpep eq "positive");
				$contro_hladrpep = "positive" if ($H =~ /$pnr/ || $pnr =~ /$H/);
			}
			
			$N_AA = $N_AA + length($pnr);	
			
			#ripristino $pnr
			$pnr =~ s/$pnr/$backup/;
# 			print "$pnr\t\t$freq_pep{$pnr}\t\t$contro_hladrpep\n";
			$pnrTABfreq = "$pnr\t\t$freq_pep{$pnr}\t\t$contro_hladrpep";
			push (@peptidixACnrTABfreq, $pnrTABfreq);
		}
		$peptidixACnrTABfreq = join('\n',@peptidixACnrTABfreq);
	}



################################################# INSERIMENTO PEPTIDI IN TABELLA MYSQL ##################################################
#	$db->query("update $Query set N_HomolPept2H_sapiens_l$min_length\_s$substitution\_m$mismatch=$NpeptidixACnr where ID='$ID'")||warn "Insertion in mysql N_HomolPept2H_sapiens_l$min_length\_s$substitution\_m$mismatch column failed for protein $ID:$!";
	
#	$db->query("update $Query set HomolPept2H_sapiens_l$min_length\_s$substitution\_m$mismatch='$peptidixACnrTABfreq' where ID='$ID'")||warn "Insertion in mysql HomolPept2H_sapiens_l$min_length\_s$substitution\_m$mismatch column failed for protein $ID:$!";

#	$db->query("update $Query set N_AA2_H_sapiens=$N_AA where ID='$ID'")||warn "Insertion in mysql column N_AA2_H_sapiens failed for protein $ID:$!";
#}


################################################ SCRITTURA DELLE DIVERSE HUMAN ENTRIES ##################################################

$H_sapiensestratte = keys(%freq_H_sapiens);
open(F,"$humandb_dir/H_sapiens.fasta");
while(<F>){
	$H_sapienstotali++ if(/^>/);
}


open (P,"+>$output_dir/$filename.matches")||die "$!";

print P "\n";
print P "##################################################################################################\n";
print P " The number of human entries mined at least once by autoimmunity.pl is $H_sapiensestratte/$H_sapienstotali\n";
print P " See below how many times each of them has been extracted                                   \n";
print P "##################################################################################################\n"; 
print P "\n";

foreach $H_sapiens (sort { $freq_H_sapiens{$b} <=> $freq_H_sapiens{$a} } keys(%freq_H_sapiens)){
	printf P ("%-95s %1d\n", $H_sapiens, $freq_H_sapiens{$H_sapiens});
}







sub contagap(){
	($stringona1_2,$z) = @_;
	$amonte = substr ($stringona1_2,0,$z);
	@gap = ($amonte =~ /(-)/g);
	$gap = @gap; 
}





$orafine = qx{date +%H:%M:%S};
print "\nautoimmunity.pl started at:\t$orainizio";
print "autoimmunity.pl ended at:\t$orafine";
