#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

################################################################################
### POD Documentation
################################################################################

=head1 NAME

    align_tophat_stats.pl  - script to calculate alignment statistics.

=head1 SYNOPSIS

    align_tophat_stats.pl  --i <path to bam list file>   
                           [--o outdir] [--v]

    parameters in [] are optional
    do NOT type the carets when specifying options

=head1 OPTIONS
    
    --i  <bam list file>   = /list of merge bam file.
   
	--b <split_bam file_list> = path to split bam file
 
   --o <output dir>       = /path/to/output directory. Optional.[PWD]

    --v                    = generate runtime messages. Optional

=head1 DESCRIPTION

The script generates summary of alignment statistics from Tophat output.

=head1 AUTHOR

Weizhong Chang,  Priti Kumari
 Bioinformatics Software Engineer I
 Institute for Genome Sciences
 University of Maryland
 Baltimore, Maryland 21201

=cut

################################################################################

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Spec;

##############################################################################
### Constants
##############################################################################

use constant FALSE => 0;
use constant TRUE  => 1;

#use constant BIN_DIR => '/usr/local/bin';


use constant PROGRAM => eval { ($0 =~ m/(\w+\.pl)$/) ? $1 : $0 };

##############################################################################
### Globals
##############################################################################

my %hCmdLineOption = ();
my $sHelpHeader = "\nThis is ".PROGRAM."\n";

GetOptions( \%hCmdLineOption,
            'outdir|o=s', 'infile|i=s', 'splitbam|b=s', 
            'verbose|v',
            'debug',
            'help',
            'man') or pod2usage(2);

## display documentation
pod2usage( -exitval => 0, -verbose => 2) if $hCmdLineOption{'man'};
pod2usage( -msg => $sHelpHeader, -exitval => 1) if $hCmdLineOption{'help'};

## make sure everything passed was peachy
check_parameters(\%hCmdLineOption);

my ($sOutDir,$prefix, $sPrefix, $sfile, $mapstats);
my ($mergebam_file, $bam_file, $mapstat_list, $mapstat_file, $f1, $path,$sf1, $sPath, $key,$pipeline1,$pipeline2);
my %bam; 
my ($cfile,$mfile,$lfile,$rfile,$pfile,$readcount,@arr,@arr1,$p_paired,$left_count,$tot_reads,$p_mapped,$fout,$out_all);
my $right_count=0;
my $bDebug   = (defined $hCmdLineOption{'debug'}) ? TRUE : FALSE;
my $bVerbose = (defined $hCmdLineOption{'verbose'}) ? TRUE : FALSE;

################################################################################
### Main
################################################################################

($bDebug || $bVerbose) ? 
	print STDERR "\nProcessing $hCmdLineOption{'infile'} ...\n" : ();

$sOutDir = File::Spec->curdir();
if (defined $hCmdLineOption{'outdir'}) {
    $sOutDir = $hCmdLineOption{'outdir'};

    if (! -e $sOutDir) {
        mkdir($hCmdLineOption{'outdir'}) ||
            die "ERROR! Cannot create output directory\n";
    }
    elsif (! -d $hCmdLineOption{'outdir'}) {
            die "ERROR! $hCmdLineOption{'outdir'} is not a directory\n";
    }
}

$sOutDir = File::Spec->canonpath($sOutDir);


####Obtaining mapstats list file from the tophat list file. 
$mapstat_list = $hCmdLineOption{'infile'};
@arr = split (/\//,$hCmdLineOption{'infile'});
$f1 = scalar @arr;
$pipeline1 = $arr[$f1-2];
$pipeline1 =~ s/_list//;
$pipeline2 = $pipeline1."_alignment_stats";
$mapstat_list=~s/bam_merge\/$pipeline1.*_list\/bam\_merge\.file\.list/samtools\_alignment\_stats\/$pipeline2\/samtools_alignment_stats\.mapstats\.list/;

open ($mergebam_file, "<$hCmdLineOption{'infile'}") or die "Error! Cannot open the tophat list file";
open ($mapstat_file, "<$mapstat_list") or die "Error! Cannot open the mapstat list file";

while (<$mergebam_file>) {
    chomp($_);
    ($f1,$path,$prefix) = File::Spec->splitpath($_);
    @arr = split(/\./,$prefix);
    $prefix = $arr[0];

    if (!(exists $bam{$prefix})) {

		open ($bam_file, "<$hCmdLineOption{'splitbam'}") or die "Error! Cannot open split bam list file";
		while(<$bam_file>){
			chomp;
#			print $_, "\n";
			if ( $_ =~ /$prefix/ ){
				($sf1,$sPath,$sPrefix) = File::Spec->splitpath($_);

				if ( -e $sPath."prep_reads.info" ) {
					$sfile = $sPath."prep_reads.info";
					open ($pfile, "<$sfile" ) or die "Error! Cannot open prep info file";
					while(<$pfile>) {
						chomp ($_);
						if ($_ =~m/left_reads_in/) {
							@arr= split (/\=/,$_);

							$bam{$prefix}{"prep_left"} += $arr[1];
						  
						}
						if ($_ =~m/right_reads_in/) {
							@arr= split (/\=/,$_);
							$bam{$prefix}{"prep_right"} += $arr[1];
		
						}
					}
					close $pfile;
 
				}

				if ( -e $sPath."left_kept_reads.info" ) {
					$sfile = $sPath.'left_kept_reads.info';
					open ($lfile, "<$sfile") or die "Error! Cannot open left info file";
					while(<$lfile>) {
						chomp ($_);

						if ($_ =~m/reads_in/) {
							@arr= split (/\=/,$_);
							$bam{$prefix}{"left"} += $arr[1];
							last;
						}
					} 
					close $lfile;
				}

	  
				if ( -e $sPath."right_kept_reads.info" ) {
					$sfile =  $sPath."right_kept_reads.info";
					open ($rfile, "< $sfile") or die "Error! Cannot open right info file";
					while(<$rfile>) { 
						chomp ($_);
						if ($_ =~m/reads_in/) {
							@arr= split (/\=/,$_);
							$bam{$prefix}{"right"} += $arr[1];
							last;
						}
					}
					close $rfile;			
				}
			}


		}
		close $bam_file;

	}
}

while (<$mapstat_file>) {

	chomp($_);
    $mapstats = $_;
#	print $_, "\n";
	($f1,$path,$prefix) = File::Spec->splitpath($mapstats);
	@arr = split(/\./,$prefix);
	$prefix = $arr[0];
#	print $prefix, "\n";
	if (exists $bam{$prefix}) {
			open ($mfile, "<$_") or die "Error! Cannot open mapstats file";
			while (<$mfile>) {
				chomp ($_);
				if ($_ =~ m/properly paired \((.*)%:/) {
					$p_paired = $1;
#					print $p_paired, "\n";
#					$p_paired =(split(/[\(:]/,$_))[1] ; 
#					$p_paired =~s/%//;
					$bam{$prefix}{'p_paired'} = $p_paired;
					last;
				}	    
		   
			}
			close $mfile;
		
			$mapstats =~ s/mapstats.txt/mapped_reads.count/;
			open ($cfile, "<$mapstats") or die "Error! Cannot open count file";
			while(<$cfile>) {
				chomp($_);
				if ($_ eq 'mapped count') {
					$readcount=<$cfile>;
					chomp($readcount);
					$bam{$prefix}{'readcount'} = $readcount;
					print 	$bam{$prefix}{'readcount'}, "\n";
				 
					last;
				}
			}
			close $cfile;
		}
}

close $mergebam_file;
close $mapstat_file; 


open ($out_all, ">$sOutDir/All_Samples.txt") or die "Error Cannot open output file";

foreach $key (keys (%bam)) {

	if (exists ($bam{$key}{'readcount'})) {
		$readcount = $bam{$key}{"readcount"};


	}

	if (exists ($bam{$prefix}{'p_paired'})){  
		$p_paired = $bam{$key}{'p_paired'};								;

	}	    
		       
	if (exists ($bam{$key}{"prep_left"})) {	

		$left_count = $bam{$key}{"prep_left"};
		
	}
	
		       
	if (exists ($bam{$key}{"prep_right"})) {	

		$left_count = $bam{$key}{"prep_right"};
		
	}

	if (exists ($bam{$key}{"left"})) {	
	
		$left_count = $bam{$key}{"left"};			
		
	}
				
	if (exists ($bam{$key}{"right"})) {

		$right_count = $bam{$key}{"right"};        
               
	}
     
	
    ###Total reads..
    $tot_reads = $left_count + $right_count;
    
    ###Percent mapped..
    $p_mapped = sprintf("%.2f",eval(($readcount/$tot_reads ) * 100)); 

    ###Write to Output Directory..
    open ($fout, ">$sOutDir/$key.txt") or die "Error Cannot open output file";
 
    print $fout "\#Sample Id\t$key\n";
    print $fout "\#Total.reads\tMapped.reads\tPercent.Mapped\tPercent.Properly.Paired\n";
    print $fout  "$tot_reads\t$readcount\t$p_mapped\t$p_paired\n";

    print $out_all "\#Sample Id\t$key\n";
    print $out_all "\#Total.reads\tMapped.reads\tPercent.Mapped\tPercent.Properly.Paired\n";
    print $out_all  "$tot_reads\t$readcount\t$p_mapped\t$p_paired\n\n";
    
    close $fout;

}   
close $out_all;

################################################################################
### Subroutines
################################################################################

sub check_parameters {
    my $phOptions = shift;
    
    ## make sure input files provided
    if (! (defined $phOptions->{'infile'}) ) {
	pod2usage( -msg => $sHelpHeader, -exitval => 1);
    }

    
}


################################################################################
