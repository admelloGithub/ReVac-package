#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
=head1  NAME 

create_table_iterator_list.pl - 

=head1 SYNOPSIS

USAGE:  create_table_iterator_list

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
                          'table_list=s', 
                          'output_iter_list=s', 
                          'log|l=s',
                          'debug=s', 
                          'help|h' ) || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Ergatis::Logger::get_logger();

if (!defined($options{'table_list'})){
    $logger->logdie("table_list was not defined");
}
if (!defined($options{'output_iter_list'})){
    $logger->logdie("output_iter_list was not defined");
}

my $list = $options{'table_list'};

## get rid of all spaces
$list =~ s/\s*//;

my @iteratorconf = split(/,/, $list);

&output_lists(\@iteratorconf, $options{'output_iter_list'});


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

    my ($iteratorconf, $output) = @_;

    open FILE, "+>$output" or $logger->logdie("Can't open output file $output");
    
    print FILE '$;TABLE$;' . "\n";

    foreach my $table (@{$iteratorconf}){

	print FILE "$table\n";
    }

    close FILE;

}
