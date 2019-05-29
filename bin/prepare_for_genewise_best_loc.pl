#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /local/packages/perl-5.8.8/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

prepare_for_geneWise.dbi

=head1 SYNOPSIS

USAGE: prepare_for_geneWise.dbi
           --project|-p                   aa1 
           --bsml_list|-i                 aat_aa.bsml.list
           --bsml_file|-I                 aat_aa_bsml_file
           --work_dir|-w                  working_directory
        [  --verbose|-v
           --PADDING_LENGTH               500
           --MIN_CHAIN_SCORE              50
           --MIN_PERCENT_CHAIN_ALIGN      70
           --num_tiers|n                  2
           --JUST_PRINT_BEST_LOCATIONS
         ]


=head1 OPTIONS

B<--project,-p>
    Name of the project database 

B<--bsml_list,-i>
    BSML list file as output by an aat_aa run

B<--bsml_file,-I>
    BSML file, a sinlge output of an aat_aa run

B<--work_dir,-w>
    Working directory to post the genome sequence segments and proteins to
    be aligned using genewise

B<--verbose,-v>
    Verbose mode 

B<--PADDING_LENGTH>
    number of basepairs to extend the genome location on each end (default: 500)

B<--MIN_CHAIN_SCORE>
    minimum score for an alignment chain to be considered a candidate for genewise realignment.
    chain_score = sum (per_id * segment_length)
    default minimum score = 50

B<--MIN_PERCENT_CHAIN_ALIGN>
    mimimum percent of the matching proteins length that aligns to the genome in a single alignment chain.
    default: 70%

B<--num_tiers,-n>
    number of overlapping best location hits.
    default: 2 (only the two best hits per location is extracted)

B<--JUST_PRINT_BEST_LOCATIONS>
    The best matches are reported to stdout.  No files are written.

B<--DUMP_CHAIN_STATS>
    All AAT alignments retrieved from the database are printed to a file \$asmbl_id.chain_stats.

B<--help,-h>
    This help documentation

=head1 DESCRIPTION

The best protein match to a given genomic region is extracted along with the genomic 
sequence corresponding to that region to be realigned using the genewise program.

The algorithm for finding the best match per location is as follows:  The results of an 
AAT search of a protein database against a genome assembly should be available for querying
from the evidence table.  All alignment chains are retrieved from the database and scored 
as the sum of (per_id * length) for each alignment segment.  The chains are sorted by score
and tiled to the genome, disallowing overlap among alignment chains along the genome.  Matches
to the top and bottom strands of the genome are tiled separately.  The first tier on each
strand contains the best hits per location.  Set --num_tiers > 1 to extract multiple best hits per
location.

Each protein is retrieved from the protein database fasta file using the cdbyank utility.  Each 
genome region is extracted as a substring of the genome sequence.  Both the protein and genome
sequence substring are written as fasta files in the working directory, structured like so:

{WORKING_DIRECTORY}/{asmbl_id}/{database}.assembly.{asmbl_id}.{counter}.pep

{WORKING_DIRECTORY}/{asmbl_id}/{database}.assembly.{asmbl_id}.{counter}.fsa

for the protein and genome sequence region, respectively.  The counter is an integer incremented for
each genome location on each asmbl_id.

The header of the .fsa file is constructed like so:

{database}.assembly.{asmbl_id}.{end5}.{end3}

so that, by extracting the header for this entry, the exact genome coordinates can be inferred.

NOTE: this script is modified version of the old database-dependent prepare_for_genewise script.

=head1 CONTACT

    Brett Whitty
    bwhitty@tigr.org

=cut

use strict;
use warnings;
use lib '/usr/local/devel/ANNOTATION/Euk_modules/bin';
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use DBI;
use Pod::Usage;
use CdbTools;
use Carp;
use XML::Twig;
use File::Basename;

#option processing
my ($database, $work_dir, $bsml_list, $bsml_file, $help, $verbose,
    $JUST_PRINT_BEST_LOCATIONS, $DUMP_CHAIN_STATS);

# default settings.
my $PADDING_LENGTH = 500; # sequence extended from each end of genomic location.
my $MIN_CHAIN_SCORE = 50; # min score of alignment required for possible genewise realignment.
my $num_tiers = 2; # only the two best matches per location is extracted.
my $MIN_PERCENT_CHAIN_ALIGN = 70; #minimum percent of the protein sequence found to align to the genome

&GetOptions (
             'project|p=s' => \$database,
             'work_dir|w=s' => \$work_dir,
             'bsml_list|i=s' => \$bsml_list,
             'bsml_file|I=s' => \$bsml_file,
             'help|h' => \$help,
             'verbose|v' => \$verbose,
             'PADDING_LENGTH=i' => \$PADDING_LENGTH,
             'MIN_CHAIN_SCORE=i' => \$MIN_CHAIN_SCORE,
             'MIN_PERCENT_CHAIN_ALIGN=f' => \$MIN_PERCENT_CHAIN_ALIGN,
             'num_tiers|n=i' => \$num_tiers,

             'JUST_PRINT_BEST_LOCATIONS' => \$JUST_PRINT_BEST_LOCATIONS,
             'DUMP_CHAIN_STATS' => \$DUMP_CHAIN_STATS,
             
         ) || pod2usage();

if ($help) {
    pod2usage();
}

unless ( $database && $work_dir && ($bsml_file || $bsml_list) 
    ) {
    carp "** Missing required options. **\n\n";
    pod2usage();
}

umask(0000); #write all output files permissibly

# track the protein alignment lengths:
my %prot_acc_to_prot_length;

## populate array of bsml files to process
my @bsml_files;
if ($bsml_list) {
    open (my $fh, $bsml_list) || die "cannot open '$bsml_list': $!";
    while (<$fh>) {
        chomp;
        if (! -e $_) {
            if (-e $_.".gz") {
                push(@bsml_files, $_.".gz");
            } else {
                die "specified bsml input file '$_' doesn't exist";
            }
        } else {
            push(@bsml_files, $_);
        }
    }
} elsif ($bsml_file) {
    if (! -e $bsml_file) {
        if (-e $bsml_file.".gz") {
            push(@bsml_files, $bsml_file.".gz");
        } else {
            die "specified bsml input file '$bsml_file' doesn't exist";
        }
    } else {
        push(@bsml_files, $bsml_file);
    }
}
if (! @bsml_files) {
    croak "error, do not have any bsml input files to process ";
}

# check for working directory:
if (! -d $work_dir) {
    mkdir($work_dir) || die "failed to create work_dir $work_dir because $!";
} 

## use XML::Twig to parse the BSML
my $twig = new XML::Twig(   
                            TwigHandlers => {
                                              'Sequence'           => \&seq_handler,
                                              'Seq-pair-alignment' => \&seq_pair_align_handler,
                                              'Analysis'           => \&analysis_handler,
                                            }
                        );

## stores alignment data (replaces db queries)
my $asmbl_data = {};
## stores assembly sequences
my $asmbl_seq = {};
## stores assembly sequence source filenames (for supporting masked aat searches)
my $asmbl_source = {};
## stores subject db path
my %subject_db;

## track the number of seq-pair-alignments we encounter
my $seq_pair_align_count = 0;

## flag for dealing properly with masked aat input BSML
my $analysis_is_masked_aat;
## path to unmasked assembly sequences (read from masked aat BSML)
my $unmasked_fasta_directory;

my $counter = 0;
my $bsml_count = scalar(@bsml_files);
foreach my $bsml_file(@bsml_files) {
    if ($verbose) {    
        print STDERR "processing ".++$counter."/$bsml_count '$bsml_file'\n";
    }
    my $infh;
    if ($bsml_file =~ /\.gz$/) {
        open($infh, "<:gzip", $bsml_file) 
         || confess("couldn't open gzipped input file '$bsml_file': $!");
    } else {
        open($infh, $bsml_file)
         || confess("couldn't open input file '$bsml_file': $!");
    }
    $twig->parse($infh);
    close $infh;

}

## if the analysis was masked aat searches we need to replace the hash of stored assembly sequence with the unmasked sequence
if ($analysis_is_masked_aat) {
    foreach my $asmbl_id(keys(%{$asmbl_source})){
        my ($source, $identifier) = split(":", $asmbl_source->{$asmbl_id});
       
        unless (-d $unmasked_fasta_directory) {
            die "unmasked fasta directory '$unmasked_fasta_directory' doesn't exist: $!";
        }
        
        $source = $unmasked_fasta_directory."/".$source;
        
        unless (-e $source) {
            die "unmasked assembly file '$source' doesn't exist: $!";
        }
        
        my $seq = get_sequence_by_id($identifier, $source);
    
        $asmbl_seq->{$asmbl_id} = $seq;
    }
}

my $protein_fasta_file;
if (scalar(keys(%subject_db)) > 1) {
    die "BSML files contained results from searches against more than one database";
} else {
    $protein_fasta_file = (keys(%subject_db))[0];
}

if (! -s $protein_fasta_file) {
    die "Subject database '$protein_fasta_file' doesn't appear to exist" if $seq_pair_align_count;
}

# get the pep and genome seqs for each location on each asmbl_id
foreach my $asmbl(keys(%{$asmbl_data})) {
    &prepare_asmbl_data($asmbl);
}

exit(0);

sub analysis_handler {
    my ($twig, $analysis) = @_;

    my $analysis_id = $analysis->{'att'}->{'id'};
    if ($analysis_id eq 'aat_aa_masked_analysis') {
        $analysis_is_masked_aat = 1;
    } else {
        return;
    }

    $unmasked_fasta_directory = get_attribute($analysis, 'unmasked_fasta_directory');
    
    unless ($unmasked_fasta_directory) {
        die "failed to parse 'unmasked_fasta_directory' Analysis Attribute";
    }

    $unmasked_fasta_directory =~ s/\/$//;
   
}   

sub seq_pair_align_handler {

    $seq_pair_align_count++;

    my ($twig, $seq_pair_alignment) = @_;

    my $seq_id = $seq_pair_alignment->{'att'}->{'refseq'};
   
    my $assembly_id = get_assembly_id($seq_id);
    my ($comp_db, $comp_id) = split(":",  $seq_pair_alignment->{'att'}->{'compxref'}, 2);
    
    if (! $comp_db || ! $comp_id) {
        die "failed parsing subject database and/or identifier from compxref '$seq_pair_alignment->{att}->{compxref}'";
    }
   
    $subject_db{$comp_db} = 1;
    
    foreach my $seq_pair_run($seq_pair_alignment->children('Seq-pair-run')) {
   
        my $chain_id = get_attribute($seq_pair_run, 'chain_number');
        my $pct_id = get_attribute($seq_pair_run, 'percent_identity');
        my $pct_sim = get_attribute($seq_pair_run, 'percent_similarity');
       
        my $orient = ($seq_pair_run->{'att'}->{'refcomplement'}) ? '-' : '+';
        
        push(@{$asmbl_data->{$assembly_id}}, 
                #$acc, $lend, $rend, $m_lend, $m_rend, $chainID, $per_id, $per_sim, $orient
             [
                $comp_id,
                $seq_pair_run->{'att'}->{'refpos'} + 1, 
                $seq_pair_run->{'att'}->{'runlength'} + $seq_pair_run->{'att'}->{'refpos'},
                $seq_pair_run->{'att'}->{'comppos'} + 1, 
                $seq_pair_run->{'att'}->{'comprunlength'} + $seq_pair_run->{'att'}->{'comppos'},
                $chain_id,
                $pct_id,
                $pct_sim,
                $orient,
             ]
            );
    }
}

####
sub prepare_asmbl_data {
    my ($asmbl_id) = @_;
    
    print "Processing assembly: $asmbl_id\n" if $verbose;
    
    ## retrieve best hits and genome locations
    my @best_hits = &get_best_location_hits($asmbl_id);
    
    if (@best_hits) {
        
        if (!$JUST_PRINT_BEST_LOCATIONS) {
            &prepare_genewise_inputs($asmbl_id, \@best_hits);
        }
    } 

    else {
        carp "warning, no best hits were returned for asmbl_id: $asmbl_id\n";
    }
}


#### 
sub get_best_location_hits {
    my ($asmbl_id) = @_;
    
    ## Populate data structure:
    
    my %alignments; #key on chainID
    # chain struct:
    #    m_lend, m_rend, lend, rend, score, accession, orient
    
    foreach my $row_ref(@{$asmbl_data->{$asmbl_id}}) {
       
        my ($acc, $lend, $rend, $m_lend, $m_rend, $chainID, $per_id, $per_sim, $orient) = @{$row_ref};
       
        print join("\t", @{$row_ref})."\n";
        
        unless ($chainID) {
            carp "Error, no chainID stored for @{$row_ref}\n";
            next;
        }
        
        $chainID = $acc . "," . $chainID;
        
        my $match_length = $m_rend - $m_lend + 1;
        
        if ($match_length < 50) { next;} # avoid shorties that falsely extend alignment chains.
        
        my $match_score = $match_length * $per_id/100;
        
        if (my $chain = $alignments{$chainID}) {
            # increment score
            $chain->{score} += $match_score;
            # adjust min/max coords for feature
            if ($chain->{lend} > $lend) {
                $chain->{lend} = $lend;
            }
            if ($chain->{rend} < $rend) {
                $chain->{rend} = $rend;
            }
            
            if ($chain->{m_lend} > $m_lend) {
                $chain->{m_lend} = $m_lend;
            }
            
            if ($chain->{m_rend} < $m_rend) {
                $chain->{m_rend} = $m_rend;
            }
        } else {
            # create chain entry:
            
            if ($lend == $rend) { next; } # don't store single coord feature
            
            $alignments{$chainID} = { lend => $lend,
                                      rend => $rend,
                                      m_lend => $m_lend,
                                      m_rend => $m_rend,
                                      score => $match_score,
                                      orient => $orient,
                                      accession => $acc,
                                      chainID => $chainID,
                                  };
        }
    }
    
    if (! %alignments) {
        return (); # no alignments found in database
    }
    
    if ($DUMP_CHAIN_STATS) {
        
        open (my $fh, ">$asmbl_id.chain_stats") or die $!;
        my @alignment_chains = sort {$a->{lend}<=>$b->{lend}} values %alignments;
        foreach my $chain (@alignment_chains) {
            my ($lend, $rend, $m_lend, $m_rend, $score, $orient, $accession, $chainID) = ($chain->{lend},
                                                                                          $chain->{rend},
                                                                                          $chain->{m_lend},
                                                                                          $chain->{m_rend},
                                                                                          $chain->{score},
                                                                                          $chain->{orient},
                                                                                          $chain->{accession},
                                                                                          $chain->{chainID});
            print $fh "$lend-$rend\[$orient]\t$m_lend-$m_rend\t$accession\t$chainID\t$score\n";
        }
        close $fh;
    }
        
    ## Store only the single best non-overlapping match for each strand
    my @top_strand_tiers;
    my @bottom_strand_tiers;
    
    ## add empty tiers:
    for (1..$num_tiers) {
        push (@top_strand_tiers, []);
        push (@bottom_strand_tiers, []);
    }

    
    my @sorted_chains = reverse sort {$a->{score}<=>$b->{score}} values %alignments;
    
    
    my @structs;
    foreach my $chain_struct (@sorted_chains) {
        
        # don't bother with low scoring alignments.
        my $score = $chain_struct->{score};
        if ($score < $MIN_CHAIN_SCORE) {
            next;
        }

        ## check the alignment length
        my ($m_lend, $m_rend) = ($chain_struct->{m_lend}, $chain_struct->{m_rend});
        my $align_len = $m_rend - $m_lend + 1;
        my $prot_acc = $chain_struct->{accession};
        

        my $prot_len = -1;
        eval {
            ## if the protein database (ie. AllGroup.niaa) has been updated since the last AAT searches
            ## we might not be able to recover certain entries from the database, in which case,
            ## this part will die!
            
            $prot_len = &get_protein_length_via_accession($prot_acc);
            
        };
        if ($@) {
            ## This is no longer tolerated.  We must be able to retrieve the relevant proteins from the fasta database.
            die "Error trying to extract protein and determine length for entry $prot_acc\n";
        }
        
        

        if ($align_len/$prot_len * 100 < $MIN_PERCENT_CHAIN_ALIGN) {
            next; #insufficient percent alignment to genome
        }
        
        ## Tile the hit
        if ($chain_struct->{orient} eq '+') {
            &try_add_chain($chain_struct, \@top_strand_tiers);
        } else {
            &try_add_chain($chain_struct, \@bottom_strand_tiers);
        }
    }
    

    
    ## move tiers to complete lists:
    my @top_strand_feats;
    foreach my $tier (@top_strand_tiers) {
        push (@top_strand_feats, @$tier);
    }
    my @bottom_strand_feats;
    foreach my $tier (@bottom_strand_tiers) {
        push (@bottom_strand_feats, @$tier);
    }
    
    if ($verbose || $JUST_PRINT_BEST_LOCATIONS) {
        print "All alignments: \n";
        my @all_alignments = sort {$a->{lend}<=>$b->{lend}} values %alignments;
        &dump_feats(@all_alignments);
        
        
        print "\n\nTop Strand Best Hits:\n";
        &dump_feats(@top_strand_feats);
        print "\n\nBottom Strand Best Hits:\n";
        &dump_feats(@bottom_strand_feats);
    }
    
    %alignments = (); #clear
    return (@top_strand_feats, @bottom_strand_feats);
}



sub prepare_genewise_inputs {
    my ($asmbl_id, $best_hits_aref) = @_;
    
    ## prepare sequence files
    my $data_dir = "$work_dir/$asmbl_id";
    
    if (! -d $data_dir) {
        mkdir ($data_dir) or croak "error, couldn't mkdir $data_dir "; 
    }
    
    my $genome_seq_length = length($asmbl_seq->{$asmbl_id});
    
    my $counter = 0;
    
    ## Prepare files: 
    foreach my $chain (@$best_hits_aref) {
        $counter++;
        
        my ($lend, $rend, $m_lend, $m_rend, $accession, $score, $orient) = ($chain->{lend},
                                                                            $chain->{rend},
                                                                            $chain->{m_lend},
                                                                            $chain->{m_rend},
                                                                            $chain->{accession},
                                                                            $chain->{score},
                                                                            $chain->{orient});
        
        my $fasta_entry;
        eval {
            $fasta_entry = cdbyank($accession, $protein_fasta_file);
        };

        if ($@) {
            ## no longer tolerated.  It's essential that this works.
            die "cdbyank error retrieving '$accession' from '$protein_fasta_file'";
        }
        

        $lend -= $PADDING_LENGTH;
        if ($lend <= 0) {
            $lend = 1;
        }
        $rend += $PADDING_LENGTH;
        if ($rend > $genome_seq_length) {
            $rend = $genome_seq_length;
        }
       
        my $startpos = $lend - 1;
        my $seq_len = $rend - $startpos;
       
        print STDERR "Getting substring($startpos, $seq_len) on assembly '$asmbl_id' which is length $genome_seq_length\n";
        
        my $subseq = substr($asmbl_seq->{$asmbl_id}, $startpos, $seq_len);
        
        if ($orient eq '-') {
            $subseq = &reverse_complement($subseq);
            ($lend, $rend) = ($rend, $lend);
        }

        $subseq =~ s/(.{1,60})/$1\n/g;

        chomp $subseq;
    
        # write genome entry
        open (my $fh, ">$data_dir/$database.assembly.$asmbl_id.$counter.fsa") or die $!;
        print $fh ">$database.assembly.$asmbl_id.$lend.$rend\n$subseq\n";
        close $fh;
        
        # write protein entry:
        open ($fh, ">$data_dir/$database.assembly.$asmbl_id.$counter.pep") or die $!;
        print $fh $fasta_entry;
        close $fh;
        
    }
}
    

####
sub try_add_chain {
    my ($chain_struct, $tier_list_aref) = @_;
    
    ## if chain_struct doesn't overlap any existing feature, go ahead and add it:
    
    
    my ($chain_lend, $chain_rend, $chain_acc, $chain_score) = ($chain_struct->{lend}, 
                                                               $chain_struct->{rend},
                                                               $chain_struct->{accession},
                                                               $chain_struct->{score});
    
  TIERS:
    foreach my $feature_list_aref (@$tier_list_aref) {    
        my $found_overlap_flag = 0;
      FEATURES:
        foreach my $feat (@$feature_list_aref) {
            my ($feat_lend, $feat_rend) = ($feat->{lend}, $feat->{rend});
            if ($feat_lend < $chain_rend && $feat_rend > $chain_lend) {
                # got overlap:
                # check to see if the accession is the same.  
                #Don't want the same protein aligned to the same region in different tiers
                if ($feat->{accession} eq $chain_acc) {
                    last TIERS;  
                }
                $found_overlap_flag = 1;
                last FEATURES;
            }
        }
        if (! $found_overlap_flag) {
            ## add to current tier:
            push (@$feature_list_aref, $chain_struct);
            last TIERS;
        }
    }
}

####
sub dump_feats {
    my @feats = @_;
    
    @feats = sort {$a->{lend}<=>$b->{lend}} @feats;
    
    foreach my $feat (@feats) {
        my ($lend, $rend, $m_lend, $m_rend, $score, $accession) = ($feat->{lend},
                                                                   $feat->{rend},
                                                                   $feat->{m_lend},
                                                                   $feat->{m_rend},
                                                                   $feat->{score},
                                                                   $feat->{accession});
        
        
        $score = sprintf ("%.2f", $score);
        print "$lend-$rend\t$accession\t$score\t$m_lend-$m_rend\n";
    }
}


sub reverse_complement { # from Egc_library.pm
    my($s) = @_;
    my ($rc);
    $rc = reverse ($s);
    $rc =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
     
    return($rc);
}

sub get_protein_length_via_accession {
    my ($prot_acc) = @_;
    
    if (my $length = $prot_acc_to_prot_length{$prot_acc}) {
        return ($length);
    }

    else {
        print "-retrieving protein length for $prot_acc\n" if $verbose;
        my $fasta_entry = cdbyank($prot_acc, $protein_fasta_file);
        my ($acc, $header, $seq) = linearize($fasta_entry);
        
        my $length = length($seq);
        $prot_acc_to_prot_length{$prot_acc} = $length;
        return ($length);
    }
}

## deals with sequence elements and makes sure that
## we have sequence for the assemblies
sub seq_handler {
    my ($twig, $sequence) = @_;
    
    my $seq_id;
    my $seq_data_import;
    my $seq_data;
    my $source;
    my $identifier;
    my $format;
    
    ## choose whether we want to deal with this sequence element
    ## based on the contents of its analysis link
    my @links = $sequence->children('Link');
    unless (@links) {
        return;
    }
    foreach my $link(@links) {
        if ($link->{'att'}->{'rel'} eq 'analysis') {
            ## may want to make role a flag option
            if ($link->{'att'}->{'role'} ne 'input_of') {
                return;
            }
        }
    }

    $seq_id = $sequence->{'att'}->{'id'};
   
    ## parse assembly id from seq id
    my $assembly_id = get_assembly_id($seq_id);
    
    if ($seq_data_import = $sequence->first_child('Seq-data-import')) {
        ## sequence is referenced via seq-data-import

        $source = $seq_data_import->{'att'}->{'source'};
        $identifier = $seq_data_import->{'att'}->{'identifier'};
        $format = $seq_data_import->{'att'}->{'format'};
        
        ## we only support pulling sequences from fasta files right now
        ## but get_sequence_by_id could be modified to support other formats 
        if (defined($format) && $format ne 'fasta') {
            confess("unsupported seq-data-import format '$format' found");
        }
        
        unless (-e $source) {
            confess("fasta file referenced in BSML Seq-data-import '$source' doesn't exist");
        }
        unless (defined($identifier)) {
            confess("Seq-data-import for '$seq_id' does not have a value for identifier");
        }
        
        $identifier =~ s/\s+//g;
        
        my $seq = get_sequence_by_id($identifier, $source);
        
        print STDERR "Retrieved sequence for '$identifier' from '$source' with length of '".length($seq)."'\n";
        
        if (length($seq) > 0) {
            $asmbl_seq->{$assembly_id} = $seq;
            ## a file with the same source basename should exist in UNMASKED_FASTA_DIRECTORY if input is masked aat output
            $asmbl_source->{$assembly_id} = basename($source).":".$identifier;
        } else {
            confess("Couldn't fetch sequence for '$seq_id' from Seq-data-import source '$source' using identifer '$identifier'");
        }
        
    } elsif ($seq_data = $sequence->first_child('Seq-data')) {
        ## sequence is in the BSML
        $asmbl_seq->{$assembly_id} = $seq_data->text();
        $asmbl_seq->{$assembly_id} =~ s/\s+//;
    } else {
        ## there is no Seq-data or Seq-data-import for the sequence
        confess("No sequence present in BSML sequence element for '$seq_id'");
    }
}

## parses out the assembly id from a sequence id 
sub get_assembly_id {
    my ($seq_id) = @_;

    if ($seq_id =~ /^[^\.]+\.assembly\.(\d+)/) {
        return $1;
    } elsif ($seq_id =~ /^_(\d+)$/) {
        return $1;
    } else {
        confess "Couldn't parse assembly id from sequence identifier '$seq_id'";
        return undef;
    }
}

sub get_attribute {
    my ($elt, $name) = @_;
    my $att = undef;

    for my $attribute ( $elt->children('Attribute') ) {
        if ( $attribute->{att}->{name} eq $name ) {
            $att = $attribute->{att}->{content};
            last;
        }
    }

    if (defined $att) {
        return $att;
    } else {
        die "failed to extract $name from attributes";
    }
}

## pull a sequence from a fasta file by sequence id
## where the sequence id is the header string up to
## the first whitespace char
sub get_sequence_by_id {
    my ($id, $fname) = @_;
    my $seq_id = '';
    my $sequence = '';
    open (IN, $fname) || confess("couldn't open fasta file for reading");
    TOP: while (<IN>) {
        chomp;
        if (/^>([^\s]+)/) {
            $seq_id = $1;
            if ($seq_id eq $id) {
                while (<IN>) {
                    chomp;
                    if (/^>/) {
                        last TOP;
                    } else {
                        $sequence .= $_;
                    }
                }
            }
        }
    }
    close IN;

    $sequence =~ s/\s+//g;
    
    return $sequence;
}
