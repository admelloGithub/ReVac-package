#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Basename;
use Ergatis::Logger;
umask(0000);

#-------------------------------------------------
# GLOBALS/DEFAULTS
#-------------------------------------------------
my $logger;

my %options = &parse_options();
my $output_iter = $options{'output'};
my $pipeline_name = $options{'pipeline_name'};

## If our tag list is not provided we create an empty
## iterator file so that the download_tag component 
## does not error out but does not download anything.
if ( defined($options{'tag_list'}) ) {
    my $tag_list = $options{'tag_list'};
    my $tags_to_download = &parse_tag_list($tag_list);
    print_download_tag_iterator($tags_to_download, $output_iter);
} else {
    create_empty_tag_list($output_iter);
}

#############################################################
####                  SUBROUTINES                        ####
#############################################################

#--------------------------------------------------
# print output iterator list that will be used 
# by ergatis
#--------------------------------------------------
sub print_download_tag_iterator {
    my ($tags, $output_iter) = @_;
    
    open (ITEROUT, "> $output_iter") or $logger->logdie("Could not open output iterator $output_iter for writing: $!");
    print ITEROUT "\$;I_FILE_BASE\$;\n";    ## Print the header identification line
    foreach my $tag (@$tags) {
        print ITEROUT "$tag\n";
    }
    
    close ITEROUT;
}

#--------------------------------------------------
# parse tag list input file and return all tag
# names to be downloaded.
#--------------------------------------------------
sub parse_tag_list {
    my $tag_list = shift;
    my $ret_tags = ();
    
    my @tag_list = split(/,/, $tag_list);
    foreach my $tag (@tag_list) {
        $tag = trim($tag);
        $tag = $pipeline_name . "_" . $tag;
        push (@$ret_tags, $tag);
    }
    
    return $ret_tags;
}

#--------------------------------------------------
# create an empty iterator file if no tags should
# be downloaded.
#--------------------------------------------------
sub create_empty_tag_list {
    my $output_iter = shift;

    open (ITEROUT, "> $output_iter") or $logger->logdie("Could not open output iterator $output_iter for writing: $!");
    print ITEROUT "\$;I_FILE_BASE\$;\n";    ## Print the header identification line
    close ITEROUT;
}

#--------------------------------------------------
# trim white space from start and end of string
#--------------------------------------------------
sub trim {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

#--------------------------------------------------
# parse command line arguments
#--------------------------------------------------
sub parse_options {
    my %opts = ();
    GetOptions(\%opts,
               'tag_list|i=s',
               'output|o=s',
               'pipeline_name|p=s',
               'log|l:s',
               'debug|d:s',
               'help') || pod2usage();
               
    if ( $opts{'help'} ) { pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ); }
    ## Set logger
    my $logfile = $opts{'log'} || Ergatis::Logger::get_default_logfilename();
    $logger = new Ergatis::Logger( 'LOG_FILE'   =>  $opts{'log'},
                                   'LOG_LEVEL'  =>  $opts{'debug'} );
    $logger = Ergatis::Logger::get_logger();
    
    ## We need to make sure our tag_list parameter and output parameter are defined
    defined ($opts{'output'}) || $logger->logdie("Please provide a proper output file iterator path");
    defined ($opts{'pipeline_name'}) || $logger->logdie("Please provide a valid pipeline name");

    return %opts;
}
