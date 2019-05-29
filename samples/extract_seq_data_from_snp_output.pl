#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";
use strict;
use warnings;

use File::Basename;

use Getopt::Long qw(:config no_ignore_case);
use IO::File;
use MUMmer::SnpDataType;

my $in		= new IO::File->fdopen(fileno(STDIN), "r");
my $out		= new IO::File("|sort|uniq>/dev/stdout");
my $INDEL_CHAR	= '.';

&parse_options;
&extract_data;

sub print_usage
{
	my $progname = basename($0);
	die << "END";
usage: $progname [--input|-i <show-snps_output>]
	[--output|-o <output>] [--help|-h]
END
}

sub parse_options
{
	my %opts = ();
	GetOptions(\%opts, "input|i=s", "output|o=s", "help|h");
	print_usage if $opts{help};
	$in->open($opts{input}, "r")
		or die "Error reading from $opts{input}: $!"
		if $opts{input};
	$out->open("|sort|uniq>$opts{output}")
		or die "Error writing to $opts{output}: $!"
		if $opts{output};
}

sub extract_data
{
	while (my $line = <$in>) {
		chomp $line;
		if ($line =~ /^\[/) {
			MUMmer::SnpDataType::SetHeader($line);
			next;
		}
		my @tokens = split /\t/, $line;
		next if scalar(@tokens) !=
			MUMmer::SnpDataType::GetNumberOfColumns();
		my $snp = new MUMmer::SnpDataType($line);
		$out->printf("%s\t%d\t%s\t%s\n", $snp->GetQueryId,
			     $snp->GetQueryPosition,
			     $snp->GetQueryFrame > 0 ? "+" : "-",
			     $snp->GetQuerySubstitution)
			if $snp->GetQuerySubstitution ne $INDEL_CHAR;
		$out->printf("%s\t%d\t%s\t%s\n", $snp->GetSubjectId,
			     $snp->GetSubjectPosition,
			     $snp->GetSubjectFrame > 0 ? "+" : "-",
			     $snp->GetSubjectSubstitution)
			if $snp->GetSubjectSubstitution ne $INDEL_CHAR;
	}
}
