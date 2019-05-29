#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
if 0; # not running under some shell

#BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
#use lib (@INC,$ENV{"PERL_MOD_DIR"});
#no lib "$ENV{PERL_MOD_DIR}/i686-linux";
#no lib ".";

#########################################
# POD DOCUMENTATION			#
#########################################

=head1 NAME

    rotate_molecule.pl  -  rotates a molecule and maps all features and evidence to new coordinate system given a coordinate as origin of replication

=head1 USAGE

    rotate_molecule.pl -d <database> -u <username> -p <password> -h <host> -a <asmbl> -s <coordinate> [ -r ]

=head1 OPTIONS

    -d <database>      database to use (Required)

    -u <username>      database username (Required)

    -p <password>      database password (Required)

    -t <host>	       database host (Required)

    -a <asmbl>         assembly to rotate (Required)

    -c <coordinate>    coordinate on current molecule that will be coordinate 1
                       on new molecule. Origin of replication  (Required)

    -reverse 		reverse complement the sequence (Optional)

    -b                 debug level. Use a large number to turn on verbose debugging (Optional)

=head1  DESCRIPTION

The program rotates an assembly given a coordinate as the origin of replication. It will then map the features on the current molecule (everything in featureloc)
to the new molecule, and adjust the genome coords stored in the table. It fetches the feature coordinates from the current assembly and updates/overwrites the featureloc
table with the newly transformed coordinates. It also changes the assembly sequence accordingly in feature table.

NOTE:  When choosing a coordinate for the --coordinate option, do note that the 5' coordinates in Manatee are 1 position greater than the fmin coordinate position in the Chado database.  The 3' and fmax positions are both the same in Manatee and Chado.  So when choosing your value, choose the Manatee coordinate.

=head1  INPUT

    None.

=head1  OUTPUT

    None

=head1  CONTACT

    Sonia Agrawal
    sagrawal@som.umaryland.edu

=begin comment

    keywords: rotate, reverse complement, molecule transformation

=end comment

=cut

use strict;
use warnings;
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use POSIX;
BEGIN {
	use Ergatis::Logger;
}

#################################################
# GLOBAL VARIABLES				#
#################################################
my %options;
my ($results, $dbh,$logger, $logfile, $asmbl_len, $asmbl_id);
my $dbtype = "mysql";

#################################################
# MAIN PROGRAM					#
#################################################
$results = GetOptions(\%options,
		'database|d=s',
		'username|u=s',
		'password|p=s',
		'host|t=s',
		'asmbl|a=s',
		'coord|c=i',
		'reverse',
		'log|l=s',
		'debug|b=s',
		'help|h') || pod2usage();

## Display documentation
if( $options{'help'} ){
	pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## Getting the log file
$logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
$logger = new Ergatis::Logger('LOG_FILE'=>$logfile,'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();
## Opening logfile for writing
open(LOGFILE, "> $logfile") or $logger->logdie("Could not open $logfile file for writing");

## Make sure everything passed was correct
&check_parameters();

if($options{'reverse'}) {
	print "Reversing ..\n";
}

## Database connection
$dbh = DBI->connect("dbi:$dbtype:host=$options{'host'};database=$options{'database'}", $options{'username'}, $options{'password'});
$logger->logdie("Database connection could not be established") if (! defined $dbh);

# Decrement the coordinate by 1 so we can do calculations correct in the 'transform' subroutine
my $offset_coord = $options{'coord'} - 1;

## Handle molecule and change the sequence to rotate it
($asmbl_len, $asmbl_id) = &build_and_load_new_molecule();

## Update all feature coordinates according to rotated molecule
&update_feature_coords();

## Closing log file
close(LOGFILE);

## Disconnect DB connection
$dbh->disconnect;

exit(0);


#################################################
# SUBROUTINES					#
#################################################

## Subroutine to change the assembly sequence based on the supplied coordinate
sub build_and_load_new_molecule {
# Get the assembly sequence from database
	my ($seq, $feat_id) = $dbh->selectrow_array("SELECT residues,feature_id FROM feature WHERE name='$options{'asmbl'}'");
	my $seqlen = length($seq);

	my $coord = $offset_coord;
	#my $coord = $options{'coord'};
# Coordinate passed should be less than the length of the assembly sequence
	if($coord > $seqlen) {
		$logger->logdie("ERROR :: Invalid coordinate supplied for rotation. Coordinate passed exceeds the length of the molecule.");
	}
	# First part is coord 0 -> input coord -1
	my $first_part = substr($seq, 0, $coord);
	# Second part is the rest
	my $second_part = substr($seq, $coord);

	# print "SEQ: " . substr($seq, 0, 10), "\n";
	# print "FIRST: " . substr($first_part, 0, 10), "\n";
	# print "SECOND: " . substr($second_part, 0, 10), "\n";
# Rotate the sequence
	my $rot_seq = $second_part . $first_part;
	if ($options{'reverse'}) {
		$rot_seq = reverse($rot_seq);
		$rot_seq =~ tr/ACGTacgt/TGCAtgca/;
    	}
	# print "SEQLEN: " . $seqlen . "\n";
	# print "ROTSEQ: " . length($rot_seq). "\n";

# Load rotated sequence into DB
	my $sth;
	$sth = $dbh->prepare("UPDATE feature SET residues = '$rot_seq' WHERE feature_id = $feat_id") or $logger->logdie($sth->errstr);
	my $rv = $sth->execute or $logger->logdie($sth->errstr);
	return ($seqlen,$feat_id);
}

###########################################################################################################################################################

## Subroutine to update the feature coordinates on the assembly molecule
sub update_feature_coords {
	my ($new_strand, $temp);
	my $featloc_hash = $dbh->selectall_hashref("SELECT featureloc_id, fmin, fmax, strand FROM featureloc WHERE srcfeature_id = $asmbl_id", 'featureloc_id');
	while (my ($flid, $ref) = each %$featloc_hash) {

		my $new_fmin = &transform($ref->{'fmin'});
		my $new_fmax = &transform($ref->{'fmax'});
		my $new_strand = $ref->{'strand'};
		my $sth;
		if($options{'reverse'}) {
			if($new_fmax < $new_fmin) {
				$temp = $new_fmin;
				$new_fmin = $new_fmax;
				$new_fmax = $temp;
			}
			if($ref->{'strand'} == -1) {
				$new_strand = 1
			} elsif($ref->{'strand'} == 1) {
				$new_strand = -1;
			}
		}

		$sth = $dbh->prepare("UPDATE featureloc SET fmin = $new_fmin, fmax = $new_fmax, strand = $new_strand WHERE featureloc_id = $flid") or $logger->logdie($sth->errstr);
		my $rv = $sth->execute or $logger->logdie($sth->errstr);
	}
}

#############################################################################################################################################################

## Subroutine to change the coordinate
sub transform {
	my $x = shift;
	my $r;
	if ($options{'reverse'}) {
		$r = ($x - $options{'coord'});
		$r = $r <= 0 ? ($r + $asmbl_len) : $r;
		$r = $asmbl_len - $r;
	} else {
		$r = $x - $offset_coord;
		$r = $r < 0 ? ($r + $asmbl_len) : $r;
	}
	return($r);
}

#################################################################################################################################################################

## Subroutine to check the supplied paramaters are correct
sub check_parameters {
	my @req = qw(database username password host asmbl coord);
	foreach my $opt ( @req ) {
		if( !$options{$opt} ) {
			$logger->logdie("Required paramater $opt not passed !");
			die("ERROR :: Required paramater $opt not passed or the parameter datatype is not correct.\n");
		}
	}
	$logger->logdie("ERROR :: Incorrect coordinate passed. Negative coordinate not allowed.") if ($options{'coord'} < 0);
	print LOGFILE "Assembly $options{asmbl} rotated at origin of replication as $options{coord}\n";
}

__END__
