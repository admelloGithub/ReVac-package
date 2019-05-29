#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

# This script is used to retrieve quota information from the
# Panasas storage given directory paths as argumetns:
#
#  Example: get_quota.pl <dir1>
#           get_quota.pl <dir1> <dir2> ... <dirN>

# Author: Victor Felix <vfelix@som.umaryland.edu>

use strict;
use IGS::Storage::Panasas::Quota;

foreach my $dir (@ARGV) {
    my $stats = get_volume_stats($dir);
    print "----------\n";
    if (defined $stats->{'Status'} && defined $stats->{'Name'}) {
        print "Name:       " . $stats->{'Name'} . "\n";
        print "Space Used: " . $stats->{'Space Used'} . "\n";
        print "Soft Quota: " . $stats->{'Soft Quota'} . "\n";
        print "Hard Quota: " . $stats->{'Hard Quota'} . "\n";
        print "Status:     " . $stats->{'Status'} . "\n";
    } else {
        print "$dir doesn't appear to be under a Panasas volume.\n";
    }
    print "----------\n";
}

