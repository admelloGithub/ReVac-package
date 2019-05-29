#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

CogBsmlLoader.pl  -  Preprocess data stored in BSML pairwise alignment documents into BTAB
structure for COG analysis using best_hits.pl. 

=head1 SYNOPSIS

USAGE:  CogBsmlLoader.pl -m BsmlGeneModelDirectory -b BsmlPairwiseAlignmentDirectory -o OutputFile

=head1 OPTIONS

=over 4

=item *

B<--bsmlModelDir, -m>   [REQUIRED] Dir containing the BSML gene/sequence encodings referenced in the search directory

=item *

B<--bsmlSearchDir, -b>  [REQUIRED] Dir containing the BSML search encodings of pairwise alignments (all_vs_all, blastp)

=item *

B<--output, -o>        [REQUIRED] output BTAB file

=item *

B<--help,-h> This help message

=back

=head1   DESCRIPTION

CogBsmlLoader.pl is designed to preprocess the data contained in a BSML pairwise alignment search 
for COGS analysis. Specifically it identifies the "best hit" per genome for each query gene. 
This data is packaged into the BTAB format for linkage analysis using best_hits.pl  

NOTE:  

Calling the script name with NO flags/options or --help will display the syntax requirement.

=cut

# Preprocess data stored in BSML pairwise alignment documents into BTAB
# structure for COG analysis.

############
# Arguments:
#
# bsmlSearchDir - directory containing BSML pairwise sequence encodings
# bsmlModelDir - directory containing the BSML gene model documents for the search data
# outfile - btab output file
#
#

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Ergatis::Logger;
use XML::Parser;

my %options = ();
my $results = GetOptions( \%options, 
              'bsmlSearchList|b=s', 
              'bsmlModelList|m=s', 
              'bsmlJaccardList|j=s', 
              'outfile|o=s', 
              'pvalcut|p=s', 
              'coverageCutoff|c=s',
              'identityCutoff|i=s',
              'similarityCutoff|s=s',
              'log|l=s',
              'debug=s',
              'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Ergatis::Logger::get_logger();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDERR} );
}


&check_parameters(\%options);

#MAIN HERE

# associative array to translate cds identifiers to polypeptide ids.
my $cds2Prot = {};

# If Jaccard data has been specified, use it for equivalence class filtering

my $jaccardClusterHash = {};  #Associate Jaccard clusters by sequence id
my $jaccardRepSeqHash = {};   #Associate a representative sequence for each cluster
my $jaccardClusterCount = 0;


$options{'bsmlJaccardList'} =~ s/\'//g;
$options{'bsmlJaccardList'} =~ s/\"//g;

if( $options{'bsmlJaccardList'} && $options{'bsmlJaccardList'} ne "" )
{
    if(-e $options{'bsmlJaccardList'}){

    my $seqCount=0;
    my $multihandlers = {'Alignment-summary'=>
                 sub {
                     $jaccardClusterCount++;
                 },
                 'Aligned-sequence'=>
                 sub {
                     my ($expat,$elt,%params) = @_;
                     my $name = $params{'name'};
                     $name =~ s/:[\d]*//;
                     # Associate a SeqId with a Jaccard Cluster ID
                     $jaccardClusterHash->{$name} = $jaccardClusterCount;
                     
                     # If this is the first sequence observed in the Jaccard Cluster
                     # identify it as the representative sequence and set it as the 
                     # representative sequence for the Jaccard Cluster ID.
                     
                     if( $seqCount == 0 )
                     {
                     $jaccardRepSeqHash->{$jaccardClusterCount} = [];
                     push @{$jaccardRepSeqHash->{$jaccardClusterCount}},$name;
                     }
                     else{
                     push @{$jaccardRepSeqHash->{$jaccardClusterCount}},$name;
                     }
                     $seqCount++;
                 }
             };
    
    my $multiparser = new XML::Parser(Handlers => 
                      {
                          Start =>
                          sub {
                              #$_[1] is the name of the element
                              if(exists $multihandlers->{$_[1]}){
                              $multihandlers->{$_[1]}(@_);
                              }
                          }
                      }
                      );
    
    open JFILE, "$options{'bsmlJaccardList'}" or $logger->logdie("Can't open file $options{'bsmlJaccardList'}");
    while(my $bsmlFile=<JFILE>){
        chomp $bsmlFile;
        $logger->debug("Parsing jaccard file $bsmlFile") if($logger->is_debug());
        if (!(-e $bsmlFile) && -e "$bsmlFile.gz") {
            $bsmlFile .= ".gz";
        }
        if(-e $bsmlFile){
            my $ifh;
            if ($bsmlFile =~ /\.(gz|gzip)$/) {
                open ($ifh, "<:gzip", $bsmlFile) || die "can't read input file $bsmlFile: $!";
            } else {
                open ($ifh, "<$bsmlFile") || die "can't read input file $bsmlFile: $!";
            }
            $multiparser->parse( $ifh );
            close $ifh;
        }
        else{
        $logger->logdie("Can't read jaccard bsml file $bsmlFile");
        }
    }
    close JFILE;

    }
    else{
    $logger->logdie("Can't read jaccard list $options{'bsmlJaccardList'}");
    }
}


# polypeptide sequence identifer to genome mapping 
my $geneGenomeMap = {};
my $genome = '';
my $instrain=0;

# loop through the documents in the model directory to create the polypeptide genome map
my $genomehandlers = {
    'Organism'=>
	sub {
	    my ($expat,$elt,%params) = @_;
	    $genome = $params{'genus'}.':'.$params{'species'}.':'.$params{'strain'};
	},
    'Strain'=>
	sub {
	    my ($expat,$elt,%params) = @_;
	    $instrain=1;
	},
    'Attribute' =>
	sub {
	    my ($expat,$elt,%params) = @_;
	    if($instrain){
		$genome .= $params{'content'};
	    }
	},  
    'Feature'=>
	sub {
	    my ($expat,$elt,%params) = @_;
	    $logger->debug("Adding $params{'id'} lookup with genome $genome");
	    $geneGenomeMap->{$params{'id'}} = $genome;
	}
};

my $genomeparser = new XML::Parser(Handlers => 
                   {
                       Start =>
			   sub {
			       #$_[1] is the name of the element
			       if(exists $genomehandlers->{$_[1]}){
				   $genomehandlers->{$_[1]}(@_);
			       }
			   },
		       End =>
			   sub {
			       if($_[1] eq 'Strain'){
				   $instrain=0;
			       }
			   }
                   }
				   );

foreach my $bsmlFile (@{&get_list_from_file($options{'bsmlModelList'})}){
    
    $logger->debug("Parsing genome file $bsmlFile") if($logger->is_debug());
    if (!(-e $bsmlFile) && -e "$bsmlFile.gz") {
        $bsmlFile .= ".gz";
    }
    my $ifh;
    if ($bsmlFile =~ /\.(gz|gzip)$/) {
        open ($ifh, "<:gzip", $bsmlFile) || die "can't read input file $bsmlFile: $!";
    } else {
        open ($ifh, "<$bsmlFile") || die "can't read input file $bsmlFile: $!";
    }
    $genomeparser->parse( $ifh );
    close $ifh;
    $genome = undef;
}

#####################################

# structure for building the COGS input. For each query gene, the COGS analysis expects
# the single best scoring alignment for each reference genome. In BSML terms, COGS expects the
# highest scoring seqpair for each refseq compseq combination where all compseqs 
# are contained in the same genome. 


#  Genome A           Genome B            Genome C
#  compseqA1          compseqB1           compseqC1
#  compseqA2          compseqB2           compseqC2
#  compseqA3                              compseqC3


# If the above represent the sets of reference genes by genome. The following would 
# correspond to the expected output if A1, B2, and C1 were the best hits by genome. 

# refseq -> compseqA1
# refseq -> compseqB2
# refseq -> compseqC1

####################################
    
my $COGInput = {};

my %alnparams;
my $bestRunScore = 0;
my $isbestrun = 0;
my $bestSeqPairRun = undef;

my $alnhandlers = {'Seq-pair-alignment'=>
             sub {
             my ($expat,$elt,%params) = @_;
             #Process the previous alignment in file
             &process_alignment($alnparams{'compseq'},$alnparams{'compgenome'},
                        $alnparams{'refseq'},$alnparams{'refgenome'},
                        $bestRunScore,$bestSeqPairRun,
                        $COGInput,$jaccardRepSeqHash,$geneGenomeMap) if(keys %alnparams);

             my $compseq = $params{'compseq'};
             my $refseq = $params{'refseq'};
             my $compGenome = $geneGenomeMap->{$compseq};
             my $refGenome = $geneGenomeMap->{$refseq};
             $params{'compgenome'}=$compGenome;
             $params{'refgenome'}=$compGenome;
             $bestRunScore = 0;
             $isbestrun=0;
             $bestSeqPairRun = undef;
             %alnparams = ();
             $logger->debug("Parsing match between $params{'compseq'} $params{'compgenome'} $params{'refseq'} $params{'refgenome'} $isbestrun");


             # self-self alignments are not included 
             return if( $compseq eq $refseq );
             
             if( !( $compGenome )) 
             {
                 $logger->debug("$compseq: compseq not found in gene genome map. skipping.") if($logger->is_debug());
                 return;
             }

             if( !( $refGenome ) )
             {
                 $logger->debug("$refseq: compseq not found in gene genome map. skipping.") if($logger->is_debug());
                 return;
             }
             
             # alignments to the same genome are not included in COG clusters
             
             return if( $compGenome eq $refGenome );
             
             %alnparams = %params;
             },
           'Seq-pair-run'=>
               sub {
               my ($expat,$elt,%params) = @_;
               if(keys %alnparams){
                   my $runscore = $params{'runscore'};
                   my $runprob = $params{'runprob'};
                   if( defined $runscore && defined $runprob
                   && ($runscore > $bestRunScore) && ($runprob < $options{'pvalcut'}) ){
                   $logger->debug("$alnparams{'compseq'} $alnparams{'refseq'} using run with runscore $runscore $runprob. Previous bestrunscore $bestRunScore. pvalue cutoff $options{'pvalcut'}");
                   $bestRunScore = $runscore;
                   $logger->debug("bestrunscore $bestRunScore");
                   $bestSeqPairRun->{'reflength'} = $alnparams{'reflength'};
                   $bestSeqPairRun->{'method'} = $alnparams{'method'};
                   $bestSeqPairRun->{'compxref'} = $alnparams{'compxref'};
                   $bestSeqPairRun->{'refpos'} = $params{'refpos'};
                   $bestSeqPairRun->{'runlength'} = $params{'runlength'};
                   $bestSeqPairRun->{'comppos'} = $params{'comppos'};
                   $bestSeqPairRun->{'comprunlength'} = $params{'comprunlength'};
                   $bestSeqPairRun->{'runscore'} = $params{'runscore'};
                   $bestSeqPairRun->{'runprob'} = $runprob;
                   $isbestrun=1;
                   }
               }
               },
           'Attribute'=>
               sub {
               my ($expat,$elt,%params) = @_;
               my $index = scalar(@{$expat->{'Context'}}) - 1;
               if($isbestrun && $expat->{'Context'}->[$index] eq 'Seq-pair-run'){
                   $logger->debug("Dumping parameters for best run $bestRunScore");
                   if($params{'name'} eq 'p_value'){
                   $bestSeqPairRun->{'p_value'} = $params{'content'};
                   }
                   elsif($params{'name'} eq 'percent_identity'){
                   $bestSeqPairRun->{'percent_identity'} = $params{'content'};
                   }
                   elsif($params{'name'} eq 'percent_similarity'){
                   $bestSeqPairRun->{'percent_similarity'} = $params{'content'};
                   }
                   elsif($params{'name'} eq 'chain_number'){
                   $bestSeqPairRun->{'chain_number'} = $params{'content'};
                   }
                   elsif($params{'name'} eq 'segment_number'){
                   $bestSeqPairRun->{'segment_number'} = $params{'content'};
                   }
               }
               if($expat->{'Context'}->[$index] eq 'Seq-pair-alignment'){
                   if($params{'name'} eq 'percent_coverage_refseq'){
                       if($params{'content'} < $options{'coverageCutoff'}) {
                           %alnparams = ();
                       }
                   }
                   if($params{'name'} eq 'percent_identity'){
                       if($params{'content'} < $options{'identityCutoff'}) {
                           %alnparams = ();
                       }
                   }
                   if($params{'name'} eq 'percent_similarity'){
                       if($params{'content'} < $options{'similarityCutoff'}) {
                           %alnparams = ();
                       }
                   }
               }
           }
           };

my $alnparser = new XML::Parser(Handlers => 
                   {
                       Start =>
                       sub {
                        #$_[1] is the name of the element
                           if(exists $alnhandlers->{$_[1]}){
                           
                           $alnhandlers->{$_[1]}(@_);
                        }
                    }
                                }
                );


open( OUTFILE, ">$options{'outfile'}" ) or $logger->logdie("Can't open file $options{'outfile'}");

foreach my $bsmlFile (@{&get_list_from_file($options{'bsmlSearchList'})}){
    
    if (!(-e $bsmlFile) && -e "$bsmlFile.gz") {
        $bsmlFile .= ".gz";
    }
    
    # builds the COGS input data structure

    $logger->debug("Parsing alignment file $bsmlFile") if($logger->is_debug());
    
    %alnparams = ();
    $bestRunScore = 0;
    $isbestrun = 0;
    $bestSeqPairRun = undef;

    my $ifh;
    if ($bsmlFile =~ /\.(gz|gzip)$/) {
        open ($ifh, "<:gzip", $bsmlFile) || die "can't read input file $bsmlFile: $!";
    } else {
        open ($ifh, "<$bsmlFile") || die "can't read input file $bsmlFile: $!";
    }
    $alnparser->parse( $ifh );
    close $ifh;
    
    #Process the last alignment in file
    &process_alignment($alnparams{'compseq'},$alnparams{'compgenome'},
               $alnparams{'refseq'},$alnparams{'refgenome'},
               $bestRunScore,$bestSeqPairRun,
               $COGInput,$jaccardRepSeqHash,$geneGenomeMap) if(keys %alnparams);;
    
    # print the results

    foreach my $k1 ( keys( %{$COGInput} ) )
    {
    foreach my $k2 (keys( %{$COGInput->{$k1}}))
    {
        my $member = $COGInput->{$k1}->{$k2}->[0];
        if(exists $jaccardClusterHash->{$member}){
        $COGInput->{$k1}->{$k2}->[21] = join(',',@{$jaccardRepSeqHash->{$jaccardClusterHash->{$member}}});
        }
        print OUTFILE join("\t", @{$COGInput->{$k1}->{$k2}});
        print OUTFILE "\n";
    }
    }

    $COGInput = {};
}


sub process_alignment{
    my($compseq,$compGenome,$refseq,$refGenome,$bestRunScore,$bestSeqPairRun,$COGInput,$jaccardRepSeqHash,$geneGenomeMap) = @_;
    # 
    if( ! defined $bestSeqPairRun ){
    $logger->warn("Best run not defined");
    return;
    }
    else{
# If compseq (or refseq) is defined in a Jaccard equivalence class identify the class by
# its reference sequence. 
    
    if( defined( my $jId = $jaccardClusterHash->{$compseq} ) )
    {
        $logger->debug("Found jaccard cluster $jId for id $compseq. Using $jaccardRepSeqHash->{$jId}->[0] as cluster representative");
        $compseq = $jaccardRepSeqHash->{$jId}->[0];
    }
    
    if( defined( my $jId = $jaccardClusterHash->{$refseq} ) )
    {
        $logger->debug("Found jaccard cluster $jId for id $refseq. Using $jaccardRepSeqHash->{$jId}->[0] as cluster representative");
        $refseq = $jaccardRepSeqHash->{$jId}->[0];
    }
    
    my $lref = [];
    
    $lref->[0] = $refseq;  #query name
    $lref->[1] = '';       #date
    $lref->[2] = $bestSeqPairRun->{ 'reflength' }; #query length
    $lref->[3] = $bestSeqPairRun->{ 'method' }; #program
    $lref->[4] = $bestSeqPairRun->{ 'compxref' };
    $lref->[5] = $compseq;
    $lref->[6] = $bestSeqPairRun->{ 'refpos' };
    $lref->[7] = $bestSeqPairRun->{ 'refpos' } + $bestSeqPairRun->{'runlength'};
    $lref->[8] = $bestSeqPairRun->{ 'comppos' };
    $lref->[9] = $bestSeqPairRun->{ 'comppos' } + $bestSeqPairRun->{'comprunlength'};
    $lref->[10] = $bestSeqPairRun->{'percent_identity'};
    $lref->[11] = $bestSeqPairRun->{'percent_similarity'};
    $lref->[12] = $bestSeqPairRun->{ 'runscore' };
    $lref->[13] = $bestSeqPairRun->{'chain_number'};
    $lref->[14] = $bestSeqPairRun->{'segment_number'};
    $lref->[15] = '';
    $lref->[16] = '';
    $lref->[17] = '';
    $lref->[18] = $bestSeqPairRun->{'comprunlength'};
    $lref->[19] = $bestSeqPairRun->{'runprob' };
    $lref->[20] = $bestSeqPairRun->{'p_value'};
    
    if($refseq && $compseq){
        if(  $COGInput->{$refseq}->{$compGenome} )
        {
        if(  $COGInput->{$refseq}->{$compGenome}->[12] < $bestRunScore )
        {
            $logger->debug("$refseq match to $compGenome with score $bestRunScore is highest scoring match.  Previous high score is $COGInput->{$refseq}->{$compGenome}->[12]");
            $COGInput->{$refseq}->{$compGenome} = $lref;
        }
        }
        else
        {
        $logger->debug("$refseq match to $compGenome is first match found.");
        $COGInput->{$refseq}->{$compGenome} = $lref;
        }
    }
    }
}


sub get_list_from_file{
    my($file) = @_;
    my @lines;
    open( FH, $file ) or $logger->logdie("Could not open $file");
    while( my $line = <FH> ){
    chomp($line);
    push @lines, split(',',$line) if($line =~ /\S+/);
    }
    return \@lines;
}

sub check_parameters{
    my ($options) = @_;
    
    if(0){
    pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});    
    }
}

# Returns only those BLAST HSPs that contributed to the best Sum(P) value
# reported for each subject/query sequence pair.
#
sub filterBlastpHsps {
    my($links) = @_;
    my $linksByQuery = &groupByMulti($links, ['from_tag', 'to_tag']);
    my $result = [];

    # Aggregate all HSPs with the same subject and query sequence
    foreach my $queryId (keys %$linksByQuery) {
	my $linksBySubject = $linksByQuery->{$queryId};

	foreach my $subjId (keys %$linksBySubject) {
	    my $slinks = $linksBySubject->{$subjId};
	    my @sortedLinks = sort { $a->{'p_value'} <=> $b->{'p_value'} } @$slinks;
	    # heuristic - assume that all HSPs with the same Sum(P) as the best are contributing to that Sum(p) score
	    my $bestScore = $sortedLinks[0]->{'p_value'};
	    
	    foreach my $sl (@sortedLinks) {
		last if ($sl->{'p_value'} > $bestScore);
		push(@$result, $sl);
	    }
	}
    }
    return $result;
}

# Generalized version of groupBy 
sub groupByMulti {
    my($arrayref, $keyFields) = @_;
    my $nKeys = scalar(@$keyFields);
    my $groups = {};

    foreach my $a (@$arrayref) {
	my @keyValues = map { $a->{$_} } @$keyFields;
	my $hash = $groups;

	for (my $i = 0;$i < $nKeys;++$i) {
	    my $kv = $keyValues[$i];

	    if ($i < ($nKeys-1)) {
		$hash->{$kv} = {} if (!defined($hash->{$kv}));
		$hash = $hash->{$kv};
	    } 
	    else {
		$hash->{$kv} = [] if (!defined($hash->{$kv}));
		push(@{$hash->{$kv}}, $a);
	    }
	}
    }
    return $groups;
}

sub getAvgBlastPPctCoverage {
    my($hsps) = @_;
    $hsps = &filterBlastpHsps($hsps);
    my $sum = 0;
    my $numHsps = 0;

    # Group by query and target id
    my $hspsByQuery = &groupByMulti($hsps, ['from_tag', 'to_tag']);

    foreach my $queryId (keys %$hspsByQuery) {
	my $hspsByTarget = $hspsByQuery->{$queryId};

	foreach my $subjId (keys %$hspsByTarget) {
	    ++$numHsps;
	    my $shsps = $hspsByTarget->{$subjId};
	    my $querySeqLen = $shsps->[0]->{'from_length'};
	    my $targetSeqLen = $shsps->[0]->{'to_length'};

	    my @queryIntervals = map { {'fmin' => $_->{'from_Nterm'}, 'fmax' => $_->{'from_Cterm'}, 'strand' => 1} } @$shsps;
	    my @targetIntervals = map { {'fmin' => $_->{'to_Nterm'}, 'fmax' => $_->{'to_Cterm'}, 'strand' => 1} } @$shsps;

	    my $mergedQueryIntervals = &mergeOverlappingIntervals(\@queryIntervals);
	    my $mergedTargetIntervals = &mergeOverlappingIntervals(\@targetIntervals);

	    my $queryHitLen = 0;
	    my $targetHitLen = 0;
	    
	    map { $queryHitLen += ($_->{'fmax'} - $_->{'fmin'}); } @$mergedQueryIntervals;
	    map { $targetHitLen += ($_->{'fmax'} - $_->{'fmin'}); } @$mergedTargetIntervals;

	    $sum += $queryHitLen / $querySeqLen;
	    $sum += $targetHitLen / $targetSeqLen;
	}
    }

    return ($numHsps > 0) ? ($sum/($numHsps * 2) * 100.0) : undef;
}

# Generate a new set of intervals by merging any that overlap in the original set.
#
sub mergeOverlappingIntervals {
    my($intervals) = @_;

    # result set of intervals
    my $merged = [];

    # sort all intervals by fmin
    my @sorted = sort { $a->{'fmin'} <=> $b->{'fmin'} } @$intervals;
    
    # current interval
    my $current = undef;

    foreach my $i (@sorted) {
	# case 1: no current interval
	if (!defined($current)) {
	    $current = $i;
	} 
	# case 2: compare current interval to interval $i
	else {
	    # case 2a: no overlap
	    if ($i->{'fmin'} > $current->{'fmax'}) {   
		push(@$merged, $current);
		$current = $i;
	    } 
	    # case 2b: overlap, with $i ending to the right of $current
	    elsif ($i->{'fmax'} > $current->{'fmax'}) {
		$current->{'fmax'} = $i->{'fmax'};
	    }
	}
    }
    push(@$merged, $current) if (defined($current));
    return $merged;
}
