#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

=head1  NAME

translate_sequence.pl  - translates nt fasta or gene describing BSML into polypeptide fasta

=head1 SYNOPSIS

USAGE:  translate_sequence.pl --transeq_bin /path/to/emboss/bin/transeq -i [single_sequence.fsa|gene_document.bsml] -o /output/directory [-l /path/to/logfile.log] [--regions s1-e1,s2-e2,...] [--frame 1,2,3,F,-1,-2,-3,R,6] [--table 0] [--multifasta_output 0] [--seqs_per_fasta 100] --project project_name --id_repository /path/to/id_repository

=head1 OPTIONS

B<--input,-i>
    The input single sequence fasta file or gene describing BSML document

B<--output,-o>
    The output path.

B<--project>
    Used as the project portion of the names in all features generated.  See Ergatis::IdGenerator
    module for more details on ID generation.

B<--regions,r>
    Optional regions of sequence to translate (formatted to be compatible with emboss transeq)
    [flag is ignored for BSML input]

B<--frame,-f>
    Optional frames to translate (1, 2, 3, F, -1, -2, -3, R, 6)
    [flag is ignored for BSML input]

B<--table,-t>
    Optional translation table to use
    
B<--multifasta_output, -m>
    Optional.  If value is set to 1, multiple sequences will be written to a fasta file, as opposed to a single seq per fasta.  Default is 0;
    Currently only works for BSML input sequences
    
B<--seqs_per_fasta -s>
	Optional.  Determines number of sequences written to a fasta file.  Only works if --multifasta_output option is enabled.
	Default is 100 seqs per file

B<--cleanup
    Optional.  If set to 1, will delete any cdb files created from the transeq executable run.  Default is 0

B<--debug>
    Debug level.  Use a large number to turn on verbose debugging.

B<--log,-l>
    Log file

B<--help,-h>
    This help message

=head1   DESCRIPTION

Translates an input fasta sequence to a polypeptide fasta file using emboss transeq.
Can also produce fasta formatted polypeptide sequences for genes described in a gene 
encoding BSML document.

=head1  CONTACT
        Brett Whitty
        bwhitty@tigr.org

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Ergatis::Logger;
use Ergatis::IdGenerator;
use XML::Twig;
use Data::Dumper;
use BSML::BsmlBuilder;
use Fasta::SimpleIndexer;
use File::Basename;

my $transeq_exec;
my $transeq_flags = ' -warning 1 -error 1 -fatal 1 -die 1 -trim 1';

my %coords = ();
my %id2title = ();
my @genes = ();
my @sequence_files = (); ## Any sequence files created by translate_sequence 
                         ## are captured here to facilitate cleanup if the 
                         ## --cleanup flag is present.

my $bsml_sequences;
my $bsml_fasta_ids;
my $bsml_sequences_topology;
my $sequence_children;
my $exon_locs;
my $cds_locs;
my $cds_regions;
my $exon_frame;
my $cds_frame;
my $polypeptide_ids;
my $multifasta_flag = 0;
my $SEQS_PER_FASTA = 100;

my %options = ();
my $results = GetOptions (\%options,
                          'transeq_bin=s',
                          'input|i=s',
                          'regions|r:s',
                          'frame|f:s',
                          'table|t:i',
                          'output|o=s',
                          'id_repository=s',
                          'project=s',
                          'cleanup:0',
                          'multifasta_output|m:0',
                          'seqs_per_fasta|s:100',
                          'log|l=s',
                          'debug=s',
                          'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

my $fasta_flag = 0;
                     
if ($options{'help'}) {
    pod2usage(verbose => 2);
}

if (!$options{'id_repository'}) {
    pod2usage("must provided --id_repository");
}
if (!$options{'project'} ) {
    pod2usage("project name must be provided with --project");
}
my $id_gen = Ergatis::IdGenerator->new('id_repository' => $options{'id_repository'});

if (!$options{'transeq_bin'}) {
    pod2usage("must provide path to transeq executable with --transeq_bin");
}
$transeq_exec = $options{'transeq_bin'};

unless ($options{'input'}) {
    pod2usage("fasta or bsml input must be provided with --input");
}

if ($options{'input'}) {
    unless (-e $options{'input'}) {
        pod2usage("input file '$options{input}' does not exist");
    }
    if ($options{'input'} =~ /\.fsa$/) {
        $fasta_flag = 1;
    }
    if ($fasta_flag && count_fasta_records($options{'input'}) != 1) {
        $logger->logdie("fasta input file must contain only one sequence record");
    }
}

if (!$options{'output'}) {
    pod2usage("must specify and output directory with --output");
}
$options{'output'} =~ s/\/$//;
unless (-d $options{'output'}) {
    pod2usage("specified output path '$options{output}' is not a directory");
}

if ($options{'multifasta_output'}) {
	if ($options{'multifasta_output'} == 0 || $options{'multifasta_output'} ==1) {
		$multifasta_flag = $options{'multifasta_output'};
	} else{
		die("The --multifasta_output option must either be a 0 (single fasta output) or a 1 (multifasta output).");
	}
	$SEQS_PER_FASTA = $options{'seqs_per_fasta'};
}


if ($options{'table'}) {
    ## translation table
    ## see transeq -h for list
    $transeq_flags .= " -table $options{table}";
}


if (!$fasta_flag) {
    ## if a bsml document has been provided
    ## we're going to process it and perform
    ## transeq on the gene models it encodes
    my $input_prefix = basename($options{'input'},".bsml");

    my $id_hash;
    
    if ($options{'regions'}) { 
         if ($logger->is_debug()) {
             $logger->debug("WARNING: --regions flag is incompatible with BSML document input and will be ignored");
         }
    }   
    if ($options{'frame'}) { 
        if ($logger->is_debug()) {
            $logger->debug("WARNING: --frame flag is incompatible with BSML document input and will be ignored");
        }
    }   
    
    ## scan through BSML input doc and find all assembly sequences
    ## make sure sequence is present as Seq-data-import or Seq-data
    parse_bsml_sequences($options{'input'});

    ## scan through BSML and pull out interval locs from feature tables 
    parse_bsml_interval_locs($options{'input'});
    
    ## constrain the ranges of the exons to within the cds feature
    foreach my $transcript_id(keys(%{$exon_locs})) {
        $exon_locs->{$transcript_id} = constrain_exons_by_cds(
                                        $exon_locs->{$transcript_id}, 
                                        $cds_locs->{$transcript_id}->[0],
                                        $cds_locs->{$transcript_id}->[1],
                                                             );
    }

    my $temp_in_fsa = $options{'output'}."/temp.in.fsa";
    my $temp_out_fsa = $options{'output'}."/temp.out.fsa";
    my $temp_multifasta = $options{'output'}."/temp.out.multi.fsa";
   
    foreach my $seq_id(keys(%{$sequence_children})) {
        my $doc = new BSML::BsmlBuilder();
        my $seq = '';
        my $seq_count = 0;	# Number of processed sequences so far
        my $multifasta_count = 1;	# ID of the current multifasta file we are writing out to
        
        ## retrieve the parent sequence
        $seq = get_sequence_by_id($bsml_sequences->{$seq_id}, $bsml_fasta_ids->{$seq_id});
    
        if (length($seq) eq 0) {
           	$logger->logdie("failed to fetch the sequence for '$seq_id'");
        }
        
        ## translate all children
        foreach my $transcript_id(@{$sequence_children->{$seq_id}}) {
           	my $flags = $transeq_flags;
        
            ## get fasta substring sequence
            my $subseq = get_substring_by_exons($seq, $exon_locs->{$transcript_id});
            
            if (length($subseq) eq 0) {
               	$logger->logdie("failed to fetch the subsequence for '$transcript_id'");
            }
            	
            write_seq_to_fasta($temp_in_fsa, "$transcript_id", $subseq);
            	
            ## transeq flags
            $flags .= " -sequence $temp_in_fsa";
            $flags .= " -outseq $temp_out_fsa";
            $flags .= " -frame '$cds_frame->{$transcript_id}'";
            $flags .= " -clean";

            ## execute transeq
            eval { system($transeq_exec . $flags) };

            if ( $@ ) {
               	$logger->logdie("failed running transeq: $@");
            }

            my $seq = $doc->createAndAddExtendedSequenceN(  'id' => $transcript_id,
                                                            'molecule' => 'na',
                                                            'class' => '',
                                                         );
            $doc->createAndAddSeqDataImportN (  'seq'           => $seq,
                                                'format'        => 'bsml',
                                                'source'        => $options{'input'},
                                                'id'            => '',
                                                'identifier'    => $transcript_id,
                                             );
            $seq->addBsmlLink('analysis', '#translate_sequence', 'input_of');
            my $ft = $doc->createAndAddFeatureTable($seq);

			if ($multifasta_flag) {
			    $seq_count++;
			    #Replace sequence IDs and write to multifasta files
            	$id_hash = replace_sequence_ids($transcript_id, $temp_out_fsa, $options{'output'}, $multifasta_count, $input_prefix);
        		if ($seq_count % $SEQS_PER_FASTA == 0) {$multifasta_count++;}   	 	      	
			} else {
            	## replace sequence ids in a single fasta output file
#         	 my $out_fsa = $options{'output'}."/"."$polypeptide_ids->{$transcript_id}.fsa";
            	$id_hash = replace_sequence_ids($transcript_id, $temp_out_fsa, $options{'output'});
			}
			
            ## remove temp out file
            unlink($temp_out_fsa);
    			
            foreach my $key(keys(%{$id_hash})) {
                my $feat = $doc->createAndAddFeature($ft, $key, $key, 'polypeptide');
                $doc->createAndAddBsmlAttribute($feat, 'reading_frame', $id_hash->{$key}->{'frame'});
                $feat->addBsmlLink('analysis', '#translate_sequence', 'computed_by');
                $feat->addBsmlLink('sequence', "#${key}_seq");
                my $transeq = $doc->createAndAddExtendedSequenceN(  'id' => $key."_seq",
                                                                    'molecule' => 'aa',
                                                                    'class' => 'polypeptide',
                                                                 );
                $doc->createAndAddSeqDataImportN(   'seq'           => $transeq,
                                                    'format'        => 'fasta',
                                                    'source'        => $id_hash->{$key}->{'fsa'},
                                                    'id'            => '',
                                                    'identifier'    => $key,
                                                );
            }
        }	   	
        
    	$doc->write($options{'output'}."/$input_prefix.translate_sequence.bsml");
    } 	
   
    ## remove temp in file
    if (-e $temp_in_fsa) {unlink($temp_in_fsa);}
    if (-e $temp_multifasta) {unlink ($temp_multifasta);}
    
    if ($options{'cleanup'} == 1) {
    	print STDERR "Removing .cdb files\n";
        &cleanup_files(@sequence_files);
    }
    
} else {
    ## otherwise we're just going to run
    ## transeq on the input nt sequence
    ## to generate a polypeptide fasta file

    my $input_prefix = basename($options{'input'},".fsa");
    
    my $doc = new BSML::BsmlBuilder();
    
    my $id_hash;
    
    my $temp_out_fsa = $options{'output'}."/"."temp.polypeptides.fsa";

    if ($options{'regions'}) {
        $transeq_flags .= " -regions $options{regions}";
    }
    
    if (!$options{'frame'}) {
        $options{'frame'} = '1';
    } 
    
    $transeq_flags .= " -frame '$options{frame}'";
    $transeq_flags .= " -sequence $options{input}";
    $transeq_flags .= " -outseq $temp_out_fsa";
    
    ## execute transeq
    eval { system($transeq_exec . $transeq_flags) };

    if ( $@ ) {
        $logger->logdie("failed running transeq: $@");
    }

    my $query_id = get_sequence_id($options{'input'});
    if( $options{'project'} eq 'parse' ) {
        $options{'project'} = $1 if( $query_id =~ /^([^\.]+)/);
    }
    
#   my $out_fsa = $options{'output'}."/"."$query_id.fsa";
    $id_hash = replace_sequence_ids('', $temp_out_fsa, $options{'output'});
    
    unlink($temp_out_fsa);
    
    my $seq = $doc->createAndAddExtendedSequenceN(  'id' => $query_id,
                                                    'molecule' => 'na',
                                                    'class' => '',
                                                 );
    $doc->createAndAddSeqDataImportN (  'seq'           => $seq,
                                        'format'        => 'fasta',
                                        'source'        => $options{'input'},
                                        'id'            => '',
                                        'identifier'    => $query_id,
                                    );
    $seq->addBsmlLink('analysis', '#translate_sequence', 'input_of');
    my $ft = $doc->createAndAddFeatureTable($seq);

    foreach my $key(keys(%{$id_hash})) {
        my $feat = $doc->createAndAddFeature($ft, $key, $key, 'polypeptide');
        $doc->createAndAddBsmlAttribute($feat, 'reading_frame', $id_hash->{$key}->{'frame'});
        $feat->addBsmlLink('analysis', '#translate_sequence', 'computed_by');
        $feat->addBsmlLink('sequence', "#${key}_seq");
        my $transeq = $doc->createAndAddExtendedSequenceN(  'id' => $key."_seq",
                                                            'molecule' => 'aa',
                                                            'class' => 'polypeptide',
                                                         );
        $doc->createAndAddSeqDataImportN(   'seq'           => $transeq,
                                            'format'        => 'fasta',
                                            'source'        => $id_hash->{$key}->{'fsa'},
                                            'id'            => '',
                                            'identifier'    => $key,
                                        );
    }
    $doc->write($options{'output'}."/$input_prefix.translate_sequence.bsml");
}


exit();

## replace the sequence ids in the output fasta file
##Output multifasta file name will be ./$prefix_$fastacount.fsa
sub replace_sequence_ids {
    my ($transcript_id, $old_fsa_file, $out_dir, $fastacount, $prefix) = @_;
    
    my $ids = {};
    
    my $count = count_fasta_records($old_fsa_file);

    if ($count > 1 && $transcript_id ne '' && ! defined ($fastacount) ) {
        $logger->logdie("Can't assign polypeptide id '$polypeptide_ids->{$transcript_id}' to more than one sequence in '$old_fsa_file'");
    }

    ## maximum number of possible translated frames is 6    
    $id_gen->set_pool_size('polypeptide' => $count);
    
    open (IN, $old_fsa_file) || $logger->logdie("couldn't open '$old_fsa_file' for reading");
    my $out_fh;
    while (<IN>) {
        if (/^>[^\s]+_(\d)/) {
            my $seq_id;
            ## if we've got a transcript id, and it has a corresponding polypeptide id
            ## we'll use it for the sequence id in the fasta header
            if (defined($transcript_id) && defined($polypeptide_ids->{$transcript_id})) {
                $seq_id = $polypeptide_ids->{$transcript_id};
            } else {
            ## otherwise we'll use the id generator
                $seq_id = $id_gen->next_id(
                                          'project' => $options{'project'},
                                          'type'    => 'polypeptide'
                                         );
                $ids->{$seq_id}->{'frame'} = $1;
            }

            if (defined $fastacount) {
            	# Multifasta output
             	my $out_fsa =  "$out_dir/$prefix" . "_" . "$fastacount.fsa";  
             	$ids->{$seq_id}->{'fsa'} = "$out_fsa";               	        
                open ($out_fh, ">>$out_fsa") || $logger->logdie("couldn't open '$out_fsa' for writing");
            } else {  
            	# Single fasta output
            	my $out_fsa = "$out_dir/$seq_id.fsa"; 
            	$ids->{$seq_id}->{'fsa'} = "$out_fsa";                    
            	open ($out_fh, ">$out_fsa") || $logger->logdie("couldn't open '$out_fsa' for writing");
            }
            print $out_fh ">$seq_id\n";
        } else {
            print $out_fh $_;
        }
    }
    close IN;
    close $out_fh;

    return $ids;
}

sub count_fasta_records {
    my ($fname) = @_;
    my $count = 0;

    my $cmd = "grep -c \">\" $fname";
    $count = &run_system_cmd($cmd);
    return $count;
}

sub get_sequence_id {
    my ($fname) = @_;
    my $id = '';
    open (IN, $fname) || $logger->logdie("couldn't open '$fname' for reading");
    while (<IN>) {
        chomp;
        if (/^>([^\s]+)/) {
            $id = $1;
            last;   
        }
    }
    close IN;

    return $id;
}

sub parse_bsml_sequences {
        my ($file) = @_;
        
        my $ifh;
        if (-e $file.".gz") {
            $file .= ".gz";
        } elsif (-e $file.".gzip") {
            $file .= ".gzip";
        }

        if ($file =~ /\.(gz|gzip)$/) {
            open ($ifh, "<:gzip", $file);
        } else {
            open ($ifh, "<$file");
        }
        
        my $twig = new XML::Twig(   TwigRoots => {'Sequence' => 1, 'Feature-group' => 1},
                                    TwigHandlers => {'Sequence' => \&process_sequence}
                                );
        
        $twig->parse($ifh);
        close $ifh;
}

sub process_sequence {
    my ($twig, $sequence) = @_;

    my $seq_id;
    my $class;
    my $molecule;
    my $seq_data_import;
    my $seq_data;
    my $source;
    my $identifier = '';
    my $format;
    my $topology;
    
    $seq_id = $sequence->{'att'}->{'id'};
    
    ## skip any non-DNA molecules or anything where the ID ends in '_seq' (since these are stubs
    #   that are referenced elsewhere within the document.)
    if ( $sequence->{'att'}->{'molecule'} ne 'dna' ) {
        $logger->warn("skipping sequence $seq_id because molecule != dna\n");
        return;
    } elsif ( $sequence->{'att'}->{'id'} =~ /_seq$/ ) {
        $logger->warn("skipping sequence $seq_id because it's a stub (id ends in _seq)\n");
        return;
    }
    
    if( $options{'project'} eq 'parse') {
        $options{'project'} = $1 if( $seq_id =~ /^([^\.])/);
    }
    $class = $sequence->{'att'}->{'class'};
    $molecule = $sequence->{'att'}->{'molecule'};
    $topology = $sequence->{'att'}->{'topology'};   
    
    if (!defined($class)) {
        if ($logger->is_debug()) {$logger->debug("Sequence class of '$seq_id' was not defined");}
    }
    if (!defined($molecule) || $molecule eq 'mol-not-set') {
        if ($logger->is_debug()) {$logger->debug("Molecule type of '$seq_id' was not defined");}
    }
    
    if (defined($topology)) {
        $bsml_sequences_topology->{$seq_id} = $topology;
    } else {
        if ($logger->is_debug()) {$logger->debug("Topology of molecule '$seq_id' was not defined");}
    }   
    
    if ($seq_data_import = $sequence->first_child('Seq-data-import')) {
        
        $source = $seq_data_import->{'att'}->{'source'};
        $identifier = $seq_data_import->{'att'}->{'identifier'};
        $format = $seq_data_import->{'att'}->{'format'};
        
        unless (-e $source) {
            $logger->logdie("fasta file referenced in BSML Seq-data-import '$source' doesn't exist");
        }
        if ($identifier eq '') {
            $logger->logdie("Seq-data-import for '$seq_id' does not have a value for identifier");
        }
        
        my $seq_count = count_fasta_records($source);
        $bsml_sequences->{$seq_id} = $source;
	$bsml_fasta_ids->{$seq_id} = $identifier;
        if ($seq_count == 1) {
            unless (get_sequence_id($source) eq $identifier) {
                print STDERR "ID disagreement between BSML Seq-data-import '$identifier' and fasta file '$source'";
            }
            
        } elsif ($seq_count > 1) {
            my $nt_seq = get_sequence_by_id($source, $identifier);
            
	    	unless (length($nt_seq) > 0) {
                $logger->logdie("couldn't extract sequence for '$seq_id' from Seq-data-import source '$source'");
	    	}

        } else {
            $logger->logdie("no fasta records found for BSML Seq-data-import '$source' for sequence '$seq_id'");
        }
            
        
    } elsif ($seq_data = $sequence->first_child('Seq-data')) {
        ## sequence is in the BSML
        ## so it will be written to a fasta file in the output dir
        
        my $sequence_file = $options{'output'}."/$seq_id.fsa";
    
        write_seq_to_fasta($sequence_file, $seq_id, $seq_data->text());
        push (@sequence_files, $sequence_file);
        
        $bsml_sequences->{$seq_id} = $sequence_file;
	$bsml_fasta_ids->{$seq_id} = $seq_id;
        
    } else {
        
        ## there is no Seq-data or Seq-data-import for the sequence
        $logger->logdie("No sequence present in BSML sequence element");
    }
   
    my $featTable = $sequence->first_child('Feature-tables');
    $sequence_children->{$seq_id} = [];
    if($featTable) {
        foreach my $child ($featTable->children('Feature-group')) {
            #$feature_parent_seq->{$child->{'att'}->{'group-set'}} = $seq_id;
            
            #Only include those sequences with a transcript id
            my ($trans_fgm) = $child->find_nodes('Feature-group-member[@feature-type="transcript"]');
            push(@{$sequence_children->{$seq_id}}, $child->{'att'}->{'group-set'}) if( $trans_fgm );
        }
    }
        
    $twig->purge;
}

sub parse_bsml_interval_locs {
        my ($file) = @_;
        
        my $ifh;

        if (-e $file.".gz") {
            $file .= ".gz";
        } elsif (-e $file.".gzip") {
            $file .= ".gzip";
        }

        if ($file =~ /\.(gz|gzip)$/) {
            open ($ifh, "<:gzip", $file);
        } else {
            open ($ifh, "<$file");
        }

        my $twig = new XML::Twig(
                                 twig_roots => {
                                        'Feature'            => \&process_feat,
                                        'Feature-group'      => \&process_feat_group,
                                               }
                                );
        $twig->parse($ifh);
}


## handles Feature elements
sub process_feat {
    my ($twig, $feat) = @_;
    my $id = $feat->att('id');

    if ($feat->att('class') eq 'exon') {
        my $seq_int = $feat->first_child('Interval-loc');
        my $complement = $seq_int->att('complement');
        my ($start_pos, $end_pos) = ($seq_int->att('startpos'), $seq_int->att('endpos'));
        if ($start_pos > $end_pos) {
            $logger->logdie("Feature '$id' startpos > endpos --- BSML usage error");
        }
        push @{$coords{$id}}, $start_pos, $end_pos;
        if ($complement) {
            $exon_frame->{$id} = -1;    
        } else {
            $exon_frame->{$id} = 1;
        }
    } elsif ($feat->att('class') eq 'CDS') {
        my $seq_int = $feat->first_child('Interval-loc');
        my ($start_pos, $end_pos) = ($seq_int->att('startpos'), $seq_int->att('endpos'));
        if ($start_pos > $end_pos) {
            $logger->logdie("Feature '$id' startpos > endpos --- BSML usage error");
        }
        push @{$coords{$id}}, $start_pos, $end_pos;
    }

    $twig->purge;
}

## handles Feature-group elements
sub process_feat_group {
    my ($twig, $feat_group) = @_;
    my @exon_coords = ();

    my $feat_group_id = $feat_group->{'att'}->{'group-set'};
    die("Could not get feat_group_id from feat_group") unless( defined( $feat_group_id ) );
    my $count = 0;
    my $sum = 0;
    my $polypeptide_flag = 0;
    foreach my $child ($feat_group->children('Feature-group-member')) {
        if ($child->att('feature-type') eq 'exon') {
            my $id = $child->att('featref');
            push(@exon_coords,[$coords{$id}->[0], $coords{$id}->[1]]);
            $sum += $exon_frame->{$id};
            $count++;
        } elsif ($child->att('feature-type') eq 'CDS') {
            my $id = $child->att('featref');
            $cds_locs->{$feat_group_id} = [$coords{$id}->[0], $coords{$id}->[1]];
        } elsif ($child->att('feature-type') eq 'polypeptide') {
            if ($polypeptide_flag) {
                $logger->logdie("Feature-group '$feat_group_id' contains more than one polypeptide feature");
            }
            $polypeptide_flag = 1;
            my $id = $child->att('featref');
            $polypeptide_ids->{$feat_group_id} = $id;
        }
    }
    if (scalar @exon_coords) {
        @exon_coords = sort { $$a[0] <=> $$b[0]; } @exon_coords;
        $exon_locs->{$feat_group_id} = \@exon_coords;
        if (abs($sum/$count) != 1) {
            $logger->logdie("transcript '$feat_group_id' has some exons on both strands");
        } else {
            $cds_frame->{$feat_group_id} = $sum/$count;
        }
    }
    $twig->purge;
}


## writes sequence to a fasta file
sub write_seq_to_fasta {
    my ($file, $header, $sequence) = @_;
    
    open (OUT, ">$file") || $logger->logdie("couldn't write fasta file '$file'");
        
    $sequence =~ s/\W+//g;
    $sequence =~ s/(.{1,60})/$1\n/g;
        
    print OUT ">$header\n$sequence";
    close OUT;
}

## pulls sequence out of a fasta file by id in header
sub get_sequence_by_id {
    my ($fname, $id) = @_;

    my $fasta = new Fasta::SimpleIndexer( $fname );
    my $sequence = $fasta->get_sequence($id);
    return $sequence;
}

## prepare substring sequence using the array of exon positions
sub get_substring_by_exons {
    my ($sequence, $exon_refs) = @_;
    
    my $subseq = '';
    
    my @cds_regions = ();
    foreach my $exon_ref(@{$exon_refs}) {
        my ($startpos, $endpos) = @{$exon_ref};
        
        if ($startpos > $endpos) {
            ## at this point we should have startpos < endpos
            ## if not, this should be resolved at the bsml parsing stage
            $logger->logdie("BSML contains an exon with startpos > endpos");
        }
    
        my $len = $endpos - $startpos;
        
        if ($startpos < 0) {
        ## for circular molecules
            $subseq .= substr($sequence, $startpos, abs($startpos));
            $subseq .= substr($sequence, 0, $endpos);
        } else {
        ## otherwise
            $subseq .= substr($sequence, $startpos, $len);
        }
    }

    return $subseq;
    
}

## constrains the exons regions to within boundaries
## set by the CDS feature 
sub constrain_exons_by_cds {
    my ($exon_loc_ref, $cds_start, $cds_end) = @_;
    
    my @exon_locs = ();
    
    foreach my $exon_ref(@{$exon_loc_ref}) {
        if ($exon_ref->[1] < $cds_start) {
            next;
        }
        if ($exon_ref->[0] > $cds_end) {
            next;
        }
        if ($exon_ref->[0] < $cds_start) {
            $exon_ref->[0] = $cds_start;
        }
        if ($exon_ref->[1] > $cds_end) {
            $exon_ref->[1] = $cds_end;
        }
        push(@exon_locs, $exon_ref);
    }
    
    return \@exon_locs;
}

## Cleanup any FASTA files created by translate_sequence.
sub cleanup_files {
    my @files = @_;

	&run_system_cmd("find $options{output} -name \"*.cdb\" -exec rm -f {} \\;");
	
    foreach my $file (@files) {
        unlink($file) or $logger->warn("Could not remove file $file.");
    }
}

sub run_system_cmd {
    my $cmd = shift;
    my $res = `$cmd`;
    chomp($res);
    my $success = $? >> 8;

    unless ($success == 0) {
        $logger->logdie("Command \"$cmd\" failed.");
    }   

    return $res;
}
