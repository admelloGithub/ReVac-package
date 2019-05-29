#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

blast2btab.pl - convert a raw wu-blast output file to btab.

=head1 SYNOPSIS

USAGE: blast2btab.pl 
            --input=/path/to/input_file.raw
            --output=/path/to/output_file.btab
          [ --log=/path/to/logfile
            --debug=N
          ]

=head1 OPTIONS

B<--input,-i>
    The input raw blast output from the wu-blastp suite.

B<--output,-o>
    The file to which the parsed output will be written.

B<--split,-s>
    Split output to one file per query. Writes files with names basename --output/queryname.btab

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

This script is used to parse a raw output file from the wu-blast suite of programs into
btab format.  The Bioperl package is used to perform this parse, as BPlite does not
currently have all the functionality necessary.

=head1  INPUT

The input to this sequence is defined using the --input option.  This should point
to the raw wu-blast output file.

=head1  OUTPUT

The output is defined using the --output option.  The file created is tab-delimited and
composed of the following columns.

    1   query_name
    2   date
    3   query_length
    4   algorithm
    5   database_name
    6   hit_name
    7   qry_start
    8   qry_end
    9   hit_start
    10  hit_end
    11  percent_identity
    12  percent_similarity
    13  raw_score
    14  bit_score
    15  NULL
    16  hit_description
    17  blast_frame
    18  qry_strand (Plus | Minus)
    19  hit_length
    20  e_value
    21  p_value

=head1  CONTACT

    Joshua Orvis
    jorvis@tigr.org

=cut


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use bigint;
BEGIN {
use Ergatis::Logger;
}
use Bio::SearchIO;

my %options = ();
my $results = GetOptions (\%options, 
                          'input|i=s',
                          'output|o=s',
			  'split|s=s',
                          'log|l=s',
			  'debug|d=s',
                          'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

my $ifh;

## get a filehandle on the input
if ($options{'input'} =~ /\.(gz|gzip)$/) {
    open ($ifh, "<:gzip", $options{'input'})
      || $logger->logdie("can't open input file:\n$!");
} else {
    open($ifh, "<$options{input}") || $logger->logdie("can't read the input sequence: $!");
}

my $in = new Bio::SearchIO(-format => 'blast', 
                           -fh     => $ifh);

## open the output file:
my $outputs = {};
my $currofh;
if(exists $options{'split'} && $options{'split'}){

}
else{
    open ($currofh, ">$options{output}") || $logger->logdie("can't create output file for BLAST report: $!");
}
# parse each blast record:
while( my $result = $in->next_result ) {

    # parse each hit per record.
    while( my $hit = $result->next_hit ) {

        # a hit consists of one or more HSPs
        while( my $hsp = $hit->next_hsp ) {
            my @x;
            $x[0] = $result->query_name();
            # date
            $x[2] = $result->query_length();
            $x[3] = $hsp->algorithm();
            
            ## database name will get parsed with whitespace if its path is long
            $x[4] = $result->database_name();
            $x[4] =~ s/\s//g;
            
            $x[5] = $hit->name();
            $x[6] = $hsp->start('query');
            $x[7] = $hsp->end('query');
            my $queryStrand = $hsp->strand('query');
            if ($queryStrand == -1) {
                ($x[6], $x[7]) = ($x[7], $x[6]);
            }

            $x[8] = $hsp->start('hit');
            $x[9] = $hsp->end('hit');
            my $hitStrand = $hsp->strand('hit');
            if ($hitStrand == -1) {
                ($x[8], $x[9]) = ($x[9], $x[8]);
            }

            $x[10] = sprintf ("%.1f", $hsp->percent_identity());   

            my $similarity = $hsp->frac_conserved('total') * 100; 
            $x[11] = sprintf("%.1f", $similarity);
            $x[12] = $hsp->score();
            $x[13] = $hsp->bits();
            
            my $desc = $hit->description();
            $desc =~ s/\\t/ /g;	#encountered rare case where a '\t' was inserted in description
            $x[15] = $desc;
            $x[16] = ( ($hsp->query->frame + 1) * $hsp->query->strand); #blast frame (1, 2, 3, -1, -2, -3).

            my $strandDescript = "null";
            if ($queryStrand == 1) {
                $strandDescript = "Plus";
            } elsif ($queryStrand == -1) {
                $strandDescript = "Minus";
            }

            $x[17] = $strandDescript;
            $x[18] = $hit->length();
            $x[19] = $hsp->evalue();
            $x[20] = $hsp->pvalue();

            if( !defined( $x[20] ) ) {
                $x[20] = &calculate_pvalue( $x[19] );
            }

	    my $outline = join ("\t", @x);
	    if(exists $options{'split'} && $options{'split'}){
		my $queryname=$x[0];
		if(!exists $outputs->{$queryname}){
		    my $dname = `dirname $options{'output'}`;
		    chomp $dname;
		    my $outputname = "$dname/$queryname.btab";
		    my $cofh;
		    open ($cofh, ">$outputname") || $logger->logdie("can't create output file for BLAST report: $!");
		    $outputs->{$queryname} = $cofh;
		}
		$currofh = $outputs->{$queryname};
	    }
            print $currofh "$outline\n";
        }
    }
}
#See http://www.ncbi.nlm.nih.gov/BLAST/tutorial/Altschul-1.html
sub calculate_pvalue {
    my $evalue = shift;

    my $estimate = 0.57721566490153;
    
    #my $p = 1 - (bexp( (-1*$evalue), 4 ) );
    ( 1 - ( $estimate**(-1*$evalue) ) );
    return $evalue;
    
}

sub check_parameters {
    my $options = shift;
    
    unless (-e $options{input}) {
        $logger->logdie("input option not passed or does not exist!");
        exit(1);
    }

    unless (defined $options{output}) {
        $logger->logdie("output option not passed");
        exit(1);
    }
    
    if(0){
        pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});    
    }
}

