#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
=head1  NAME 

create_iterator_list.pl - 

=head1 SYNOPSIS

USAGE:  create_iterator_list

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
use Ergatis::Logger;
use File::Basename;

umask(0000);

my %options = ();


my $results = GetOptions (\%options, 
                          'input_list=s', 
			  'output_iter_list=s',
			  'log|l=s',
			  'key=s',
                          'debug=s', 
                          'help|h' ) || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Ergatis::Logger::get_logger();

if (!defined($options{'input_list'})){
    $logger->logdie("input_list was not defined");
}

if (!defined($options{'output_iter_list'})){
    $logger->logdie("output_iter_list was not defined");
}

my $list = $options{'input_list'};

## get rid of all spaces
$list =~ s/\s*//;

my @iteratorconf = split(/,/, $list);

&output_lists(\@iteratorconf, $options{'output_iter_list'}, $options{'key'});


print "$0 program execution completed\n";
print "Log file is '$logfile'\n";
exit(0);
						     
#---------------------------------------------------------------------------------------------------------
#
#                           END OF MAIN  --  SUBROUTINES FOLLOW
#
#---------------------------------------------------------------------------------------------------------


#---------------------------------------------
# output_lists()
#
#---------------------------------------------
sub output_lists {

    my ($iteratorconf, $output, $key) = @_;

    open FILE, "+>$output" or $logger->logdie("Can't open output file $output");
    
    if (!defined($key)){
	print FILE '$;CV_ID$;' . "\n";
    }
    else {
	print FILE '$;' . $key . '$;' . "\n";
    }

    foreach my $cv_id (@{$iteratorconf}){

	print FILE "$cv_id\n";
    }

    close FILE;

}
