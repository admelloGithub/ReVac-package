#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

=head1 NAME

ps_scan2bsml.pl - Formats prosite scan output into bsml format.

=head1 SYNOPSIS

USAGE: template.pl
            --input_file=/path/to/some/transterm.raw
            --output=/path/to/transterm.bsml
            --analysis_id=ps_scan
          [ --log=/path/to/file.log
            --debug=4
            --help
          ]

=head1 OPTIONS

B<--input_file,-i>
    The input file (should be prosite scan output)

B<--output,-o>
    Where the output bsml file should be

B<--analysis_id,-a>
    Analysis id.  Should most likely by ps_scan_analysis.
    
B<--query_file_path,-g>
    Path to the query file (input fasta file) for ps_scan.

B<--gzip_output,-g>
    A non-zero value will result in compressed bsml output.  If no .gz is on the end of the bsml output name, one will
    be added.

B<--log,-l>
    In case you wanted a log file.

B<--debug,-d>
    There are no debug statements in this program.  Sorry.

B<--help,-h>
    Displays this message.

=head1  DESCRIPTION

    Reads ps_scan output file and turns it into BSML.  Generates BSML with Sequence Pair Alignment and 
    Sequence Pair Run elements with polypeptide and prosite_entry class elements.

=head1  INPUT

    The raw output of ps_scan.  That is, the default output without using the -o option
    for ps_scan.pl.

=head1  OUTPUT

    Generates a BSML document representing the ps_scan matches.

=head1  CONTACT

    Kevin Galens
    kgalens@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use BSML::BsmlBuilder;
use Ergatis::Logger;
use Ergatis::IdGenerator;
use Data::Dumper;

####### GLOBALS AND CONSTANTS ###########
my $inputFile;
my $output;
my $debug;
my $analysis_id;
my $fasta_file;
my $gzip;
my $identifier_lookup = {};
########################################

my %options = ();
my $results = GetOptions (\%options, 
                          'input_file|i=s',
                          'output|o=s',
                          'analysis_id|a=s',
                          'query_file_path|q=s',
                          'gzip_output|g=s',
                          'fasta_file|f=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || &_pod;

#Setup the logger
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

# Check the options.
&check_parameters(\%options);

my $data = &parsePs_scanData($inputFile);
my $bsml = &generateBsml($data);
print "Writing to $output\n";
$bsml->write($output,'', $gzip);

######################## SUB ROUTINES #######################################
sub parsePs_scanData {
    my $file = shift;
    my $retHash;
    my ($seq, $prosite);

    #Open the file (should be default output from ps_scan run).
    open(IN, "<$file") or 
        &_die("Unable to open ps_scan raw file input ($file) : $!");

    while(<IN>) {

        if(/^>(.*)\s:(.*)$/) {
            $seq = $1;
            $seq =~ s/\s+//g;
            $prosite = $2;

            &_die("prosite ($prosite) or seq ($seq) was not set\n") unless($prosite && $seq);
        } else {
            my @tmp = split(/[\s\t]+/);

            if($tmp[1] > $tmp[3]) {
                $tmp[0] = 1;
                ($tmp[1], $tmp[3]) = ($tmp[3], $tmp[1]);
            } else {
                $tmp[0] = 0;
            }

            my $match = { 'start'  => $tmp[1]-1,
                          'stop'   => $tmp[3]-1,
                          'strand' => $tmp[0],
                          'match'  => $tmp[4] };

            push(@{$retHash->{$seq}->{$prosite}}, $match);
        }

    }

    return $retHash;

}

sub generateBsml {
    my $data = shift;
    my $seqObj;

    ## keeps track of which sequence elements we've added (since they can't be duplicated)
    ##  key is id, value is a reference to the object
    my %seqs;
    
    my $doc = new BSML::BsmlBuilder();

    if( scalar(keys %{$data}) == 0 ) {
        foreach my $seq ( keys %{$identifier_lookup} ) {
            print "Adding $seq to doc\n";
            my $seq_elem = $doc->createAndAddSequence($seq, $seq, '', 'aa', 'polypeptide');
            $doc->createAndAddLink($seq_elem, 'analysis', '#'.$analysis_id, 'input_of');
            $doc->createAndAddSeqDataImport($seq_elem, 'fasta', $fasta_file, '', $seq);
            $doc->createAndAddBsmlAttribute( $seq_elem, 'defline', $identifier_lookup->{$seq});
        }
    }

    foreach my $seq(keys %{$data}) {
        
        #Create the seq object if we need to
        if ( ! defined $seqs{$seq} ) {
            $seqs{$seq} = $doc->createAndAddSequence($seq, $seq, '', 'aa', 'polypeptide');
        }
        
        $doc->createAndAddLink($seqs{$seq}, 'analysis', '#'.$analysis_id, 'input_of');

        if( ! exists( $identifier_lookup->{$seq} ) ) {
            die("$seq was not found in query sequence file");
        }

        $doc->createAndAddSeqDataImport($seqs{$seq}, 'fasta', $fasta_file, '', $seq);
        $doc->createAndAddBsmlAttribute( $seqs{$seq}, 'defline', $identifier_lookup->{$seq});

        foreach my $proDomain(keys %{$data->{$seq}}) {

            #Create another seq object.
            my ($proId, $title) = ($1, $2) if($proDomain =~ /(\w+)\s(.*)/);
            &_die("Could not parse id and title from prosite domain id line in raw output")
                unless($proId && $title);
                
            if ( ! defined $seqs{$proId} ) {
                $seqs{$proId} = $doc->createAndAddSequence($proId, $title,  '', 'aa', 'prosite_entry');
            }

            my %spaArgs = ( 'refseq'  => $seqs{$seq}->{'attr'}->{'id'},
                            'compseq' => $seqs{$proId}->{'attr'}->{'id'},
                            'restart' => 0,
                            'class'   => 'match' );

            my $spaObj = $doc->createAndAddSequencePairAlignment(%spaArgs);

            foreach my $match (@{$data->{$seq}->{$proDomain}}) {
                
                my %sprArgs = ( 'alignment_pair' => $spaObj,
                                'refpos'         => $match->{'start'},
                                'runlength'      => $match->{'stop'} - $match->{'start'},
                                'refcomplement'  => $match->{'strand'},
                                'comppos'        => 0,
                                'comprunlength'  => length($match->{'match'}),
                                'compcomplement' => 0,
                                );   

                my $sprObj = $doc->createAndAddSequencePairRun(%sprArgs);
                $doc->createAndAddBsmlAttribute($sprObj, 'class', 'match_part');
            }

            ## add a link to the analysis
            my $analysis_link = $doc->createAndAddLink( $spaObj, 'analysis', "#$analysis_id", 'computed_by' );


        }
        


    }
    $doc->createAndAddAnalysis( 'id' => $analysis_id,
                                'sourcename' => $output,
                                'algorithm'  => 'ps_scan',
                                'program'    => 'ps_scan');

    return $doc;
}

sub check_parameters {
    my $options = shift;

    my $error = "";

    &_pod if($options{'help'});

  
    if($options{'input_file'}) {
        $error .= "Option input_file ($options{'input_file'}) does not exist\n" unless(-e $options{'input_file'});
        $inputFile = $options{'input_file'};
    } else { 
        $error .= "Option input_file is required\n";
    }

    unless($options{'output'}) {
        $error .= "Option output is required\n";
    } else {
        $output = $options{'output'};
    }

    unless($options{'analysis_id'}) {
        $error .= "Option analysis_id is required\n";
    } else {
        $analysis_id = $options{'analysis_id'};
    }

    unless($options{'query_file_path'}) {
        $error .= "Option fasta_file is required\n";
    } else {
        $fasta_file = $options{'query_file_path'};
        open(IN, "< $fasta_file") or
            &_die("Unable to open $fasta_file ($!)");
        
        while(<IN>) {
            if(/^>(.*)/) {
                my $defline = $1;
                my $identifier = $1 if($defline =~ /^([^\s]+)/);
                $identifier_lookup->{$identifier} = $defline;
            }
        }
        close(IN);
    }
    
    if($options{'debug'}) {
        $debug = $options{'debug'};
    }

    if($options{'gzip_output'}) {
        $gzip = 1;
    } else {
        $gzip = 0;
    }

    
    
    unless($error eq "") {
        &_die($error);
    }
    
}

sub _pod {   
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

sub _die {
    my $msg = shift;
    $logger->logdie($msg);
}
