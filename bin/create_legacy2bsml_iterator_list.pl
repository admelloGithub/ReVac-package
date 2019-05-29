#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
=head1  NAME 

generate_asmbl_list.pl - Default output is a workflow iterator that
can be used to iterator over a set of asmbl_ids

=head1 SYNOPSIS

USAGE:  generate_asmbl_list

=head1 OPTIONS

=item *

B<--debug,-d> Debug level.  Use a large number to turn on verbose debugging. 

=item *

B<--log,-l> Log file

=item *

B<--help,-h> This help message

=head1   DESCRIPTION

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

use File::Basename;
use Annotation::Util2;
use Ergatis::Logger;

umask(0000);

## These are the only valid organism types/schemas
my $valid_organism_types = { 'prok'   => 1,
			     'ntprok' => 1,
			     'euk'    => 1
			 };

## These are the headers that should be written to the output file.
my $HEADER_LIST = ['$;UNIQUE_KEY$;',
		   '$;DATABASE$;',
		   '$;ASMBL_ID$;',
		   '$;SEQUENCE_TYPE$;',
		   '$;SCHEMA_TYPE$;',
		   '$;INCLUDE_GENEFINDERS$;',
		   '$;EXCLUDE_GENEFINDERS$;',
		   '$;ALT_DATABASE$;',
		   '$;ALT_SPECIES$;',
		   '$;TU_LIST_FILE$;',
		   '$;MODEL_LIST_FILE$;',
		   '$;GOPHER_ID_MAPPING_FILE$;',
		   '$;REPEAT_ID_MAPPING_FILE$;',
		   ];
		   
# set up options
my ($controlFile, $outfile, $logfile, $help, $debug, $man,
    $gopherIdMappingFile, $repeatIdMappingFile);

&GetOptions(
	    'control_file|f=s' => \$controlFile,
	    'output|o=s'      => \$outfile,
	    'log|l=s'          => \$logfile,
	    'help|h'           => \$help,
	    'debug=s'          => \$debug,
	    'man|m'            => \$man,
	    'gopher_file=s'    => \$gopherIdMappingFile,
	    'repeat_file=s'    => \$repeatIdMappingFile,
	    );


if (!defined($logfile)){
    $logfile = Ergatis::Logger::get_default_logfilename();
}

my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$debug);

$logger = Ergatis::Logger::get_logger();


&checkCommandLineArguments();


if (! Annotation::Util2::checkInputFileStatus($controlFile)){
    $logger->logdie("Detected some problem with the control ".
		    "file '$controlFile'");
}

my $iteratorconf = [];

#
# Read in the information from asmbl_file OR asmbl_list
#
&get_list_from_file($iteratorconf,
		    $controlFile,
		    $gopherIdMappingFile,
		    $repeatIdMappingFile);

#
# Output the lists
#
&output_lists($iteratorconf, $outfile);


print "$0 execution completed\n";
print "The output file is '$outfile'\n";
print "The log file is '$logfile'\n";
exit(0);
						     
#---------------------------------------------------------------------------------------------------------
#
#                           END OF MAIN  --  SUBROUTINES FOLLOW
#
#---------------------------------------------------------------------------------------------------------

sub checkCommandLineArguments {

    if ($help){
	&pod2usage( {-exitval => 1, -verbose => 2, -output => \*STDOUT} ); 
	exit(1);
    }
    
    my $fatalCtr=0;

    if (!$controlFile){
	print STDERR "control_file was not specified\n";
	$fatalCtr++;
    }
    if (!$outfile){
	print STDERR "output was not specified\n";
	$fatalCtr++;
    }


    if ($fatalCtr>0){
	die "Required command-line arguments were not specified";
    }


}



#-----------------------------------------
# get_list_from_file()
#
#-----------------------------------------
sub get_list_from_file{

    my ($iteratorconf, $controlFile, $gopherIdMappingFile, $repeatIdMappingFile) = @_;

    if (defined($gopherIdMappingFile)){
	if (! Annotation::Util2::checkInputFileStatus($gopherIdMappingFile)){
	    $logger->logdie("Detected some problem with the GOPHER ".
			    "identifier mapping file '$gopherIdMappingFile'");
	}
    } else {
	$gopherIdMappingFile = '';
    }

    if (defined($repeatIdMappingFile)){
	if (! Annotation::Util2::checkInputFileStatus($repeatIdMappingFile)){
	    $logger->logdie("Detected some problem with the repeat ".
			    "identifier mapping file '$repeatIdMappingFile'");
	}
    } else {
	$repeatIdMappingFile = '';
    }

    my $contents = Annotation::Util2::getFileContentsArrayRef($controlFile);

    if (!defined($contents)){
	$logger->logdie("Could not retrieve contents of control file ".
			"'$controlFile'");
    }

    my $orghash = &get_organism_hash($contents, $controlFile);

    foreach my $database_type (sort keys %{$orghash} ){ 

	my $database      = $orghash->{$database_type}->{'database'};
	my $organism_type = $orghash->{$database_type}->{'organism_type'};
	my $include_genefinders = $orghash->{$database_type}->{'include_genefinders'};
	my $exclude_genefinders = $orghash->{$database_type}->{'exclude_genefinders'};
	my $alt_database = $orghash->{$database_type}->{'alt_database'};
	my $alt_species  = $orghash->{$database_type}->{'alt_species'};

	foreach my $infohash ( @{$orghash->{$database_type}->{'infohash'}} ) {

	    my $sequence_type = $infohash->{'sequence_type'};
	    my $asmbl_id = $infohash->{'asmbl_id'};
	    my $tu_list = $infohash->{'tu_list'};
	    my $model_list = $infohash->{'model_list'};

	    my $uniqueId = $database . '_' . $asmbl_id;
	    
	    push(@{$iteratorconf}, [ $uniqueId, 
				     $database, 
				     $asmbl_id,
				     $sequence_type, 
				     $organism_type,
				     $include_genefinders, 
				     $exclude_genefinders, 
				     $alt_database,
				     $alt_species,
				     $tu_list,
				     $model_list,
				     $gopherIdMappingFile,
				     $repeatIdMappingFile]);
	}
    }
}

#---------------------------------------------
# output_lists()
#
#---------------------------------------------
sub output_lists {

    my ($iteratorconf, $output) = @_;

    open FILE, "+>$output" or $logger->logdie("Can't open output file $output");
    

    print FILE join("\t",@{$HEADER_LIST}) ,"\n";

    foreach my $arrayRef (@{$iteratorconf}){

	print FILE join("\t",@{$arrayRef}),"\n";

    }

    close FILE;

}


#---------------------------------------------
# get_organism_hash()
#
#---------------------------------------------
sub get_organism_hash {

    my ($contents, $file) = @_;

    my $hash = {};

    my $database_type;

    my $unique_asmbl_id_values = {};

    my $linectr=0;

    foreach my $line (@{$contents}){

	$linectr++;

	if ($line =~ /^\s*$/){
	    next; # skip blank lines
	}
	elsif ($line =~ /^\#/){
	    next; # skip comment lines
	}
	elsif ($line =~ /^\-\-/){
	    next; # skip -- lines
	}
	else{

	    if ($line =~ /^database:(\S+)\s+organism_type:(\S+)\s+include_genefinders:(\S*)\s+exclude_genefinders:(\S*)\s+alt_database:(\S*)\s+alt_species:(\S*)\s*$/){
		
		my $database      = $1;
		my $organism_type = $2;
		my $include_genefinders = $3;
		my $exclude_genefinders = $4;

		my ($alt_database, $alt_species) = &verify_alt_database_and_species($5, $6);

		$database_type = $database ."_" .$organism_type;

		if (&verify_organism_type($organism_type, $linectr)){
		    
		    ($include_genefinders, $exclude_genefinders) = &verify_and_set_genefinders($include_genefinders, $exclude_genefinders, $linectr);

		    

		    if (( exists $hash->{$database_type}) &&
			(defined($hash->{$database_type}))){
			
			$logger->warn("This database_type '$database_type' was already encountered in the control file!");
		    }
		    else {
			#
			# Encountered this organism/database for the first time while processing the control file contents
			# therefore go ahead and declare the organism's hash attributes
			#
			$hash->{$database_type} = { 'database'            => $database,
						    'asmbl_id_list'       => [],
						    'infohash'            => [],
						    'organism_type'       => $organism_type,
						    'include_genefinders' => $include_genefinders,
						    'exclude_genefinders' => $exclude_genefinders,
						    'alt_database'        => $alt_database,
						    'alt_species'         => $alt_species
						};
			
		    }
		}
	    }
	    elsif ($line =~ /^database:(\S+)\s+organism_type:(\S+)\s+include_genefinders:(\S*)\s+exclude_genefinders:(\S*)/){

		my $database      = $1;
		my $organism_type = $2;
		my $include_genefinders = $3;
		my $exclude_genefinders = $4;

		$database_type = $database ."_" .$organism_type;

		if (&verify_organism_type($organism_type, $linectr)){
		    
		    ($include_genefinders, $exclude_genefinders) = &verify_and_set_genefinders($include_genefinders, $exclude_genefinders, $linectr);

		    

		    if (( exists $hash->{$database_type}) &&
			(defined($hash->{$database_type}))){
			
			$logger->warn("This database_type '$database_type' was already encountered in the control file!");
		    }
		    else {
			#
			# Encountered this organism/database for the first time while processing the control file contents
			# therefore go ahead and declare the organism's hash attributes
			#
			$hash->{$database_type} = { 'database'            => $database,
						    'asmbl_id_list'       => [],
						    'infohash'            => [],
						    'organism_type'       => $organism_type,
						    'include_genefinders' => $include_genefinders,
						    'exclude_genefinders' => $exclude_genefinders,
						    'alt_database'        => 'none',
						    'alt_species'         => 'none'
						};
			
		    }
		}
	    }
	    elsif ($line =~ /^database:(\S+)\s*/){
		$database_type = $1;
		my $database = $database_type;
		$hash->{$database_type}->{'database'} = $database;
	    }
	    elsif ($line =~ /^\s*(\d+)\s*/){

		&store_asmbl_id($database_type, $1, $line, $linectr, $unique_asmbl_id_values, $hash, $file);
	    }
	    else {
		$logger->logdie("Could not parse line number '$linectr' - line was '$line'");
	    }
	}
    }

    return $hash;

}



#------------------------------------------------------
# store_asmbl_id()
#
#------------------------------------------------------
sub store_asmbl_id {
    
    my ($database, $asmbl_id, $line, $linectr, $unique_asmbl_id_values, $hash, $file) = @_;

    if (( exists $unique_asmbl_id_values->{$database}->{$asmbl_id}) && 
	(defined($unique_asmbl_id_values->{$database}->{$asmbl_id}))){
	
	$logger->logdie("Already processed information for asmbl_id '$asmbl_id' organism '$database'.  Please review legacy2bsml control file '$file'.");
    }
    else {

	my $sequence_type = 'none';
	
	my $tu_list_file = 'none';

	my $model_list_file = 'none';

	my @attributes = split(/\s+/, $line);
	
	
	foreach my $attribute (@attributes){

	    if ($attribute =~ /:/){

		my ($key, $value) = split(/:/, $attribute);
		
		if ($key eq 'sequence_type'){
		    $sequence_type = &verify_and_set_sequence_type($value, $linectr);
		}
		elsif ($key eq 'tu_list_file'){
		    $tu_list_file = &verify_and_set_tu_list_file($value);
		}
		elsif ($key eq 'model_list_file'){
		    $model_list_file = &verify_and_set_tu_list_file($value);
		}
		else {
		    $logger->warn("Unrecognized attribute");
		}
	    }
	}

	## Store the next valid asmbl_id in the organism hash
	my $infohash = { 'asmbl_id'      => $asmbl_id,
			 'sequence_type' => $sequence_type,
			 'tu_list'       => $tu_list_file,
			 'model_list'    => $model_list_file 
		     };
	
	push( @{$hash->{$database}->{'infohash'}}, $infohash);

	## Update the unique_asmbl_id_values hash to ensure unique assembly identifiers are processed
	$unique_asmbl_id_values->{$database}->{$asmbl_id}++;

    }
}


#-----------------------------------------------------
# verify_and_set_tu_list_file()
#
#-----------------------------------------------------
sub verify_and_set_tu_list_file {

    my ($file) = @_;

    if (-e $file){
	if (-f $file){
	    if (-r $file){
		if (!-z $file){
		    return $file;
		}
		else {
		    $logger->logdie("file '$file' had zero content");
		}
	    }
	    else {
		$logger->logdie("file '$file' does not have read permissions");
	    }
	}
	else {
	    $logger->logdie("file '$file' is not a regular file");
	}
    }
    else {
	$logger->logdie("file '$file' does not exist");
    }
}


#-----------------------------------------------------
# verify_organism_type()
#
#-----------------------------------------------------
sub verify_organism_type {

    my ($organism_type, $line) = @_;
    
    if ((exists $valid_organism_types->{$organism_type}) && 
	(defined($valid_organism_types->{$organism_type}))){
	
	return 1;
    }
    else {
	$logger->logdie("Encountered an invalid organism type '$organism_type' at line '$line'");
    }

}

#-----------------------------------------------------
# verify_and_set_sequence_type()
#
#-----------------------------------------------------
sub verify_and_set_sequence_type {

    my ($sequence_type, $linectr) = @_;
    
    
    if ((!defined($sequence_type)) or ($sequence_type =~ /^\s*$/)) {
	$sequence_type = "none";
    }
    
    #
    # Lame. That's all the checking for now.
    #
    return ($sequence_type);
    
    
}

#-----------------------------------------------------
# verify_and_set_genefinders()
#
#-----------------------------------------------------
sub verify_and_set_genefinders {

    my ($include_genefinder, $exclude_genefinder, $linectr) = @_;
    
    if ((!defined($exclude_genefinder)) || 
	($exclude_genefinder =~ /^\s*$/)){

	$exclude_genefinder = 'none';
    } 
    if ((!defined($include_genefinder)) ||
	($include_genefinder =~  /^\s*$/)){

	$include_genefinder = 'all';  
    }

    if ($include_genefinder eq 'none'){
	$exclude_genefinder = 'all';
    }
    elsif ($exclude_genefinder eq 'all'){
	$include_genefinder = 'none';
    }
        
    return ($include_genefinder, $exclude_genefinder);
    
    
}


sub verify_alt_database_and_species {
    
    my ($alt_database, $alt_species) = @_;
 
    if ((!defined($alt_database)) ||
	($alt_database =~ /^\s*$/)){
	
	$alt_database = "none";
    }

    if ((!defined($alt_species)) ||
	($alt_species =~ /^\s*$/)){
	
	$alt_species = "none";
    }

    return ($alt_database, $alt_species);
}
