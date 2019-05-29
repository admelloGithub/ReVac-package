package Split_DB_File;

use warnings;
use strict;

use DB_File;

require Exporter;

our @ISA = qw(Exporter);

our @EXPORT = @DB_File::EXPORT;

sub TIEHASH
{
	my $class = shift;
	my $this = {};
	$this->{db_size} = 500000;
	$this->{db_name} = shift || die "Database filename is required";
	$this->{dbs} = [];
	$this->{idx} = 0;
	@{$this->{opts}} = @_;
	bless $this, $class;
	$this->init_hash();
	return $this;
}

sub FETCH
{
	my $this = shift;
	my ($key) = @_;
	my $db = $this->find_db_by_key($key);
	return $db ? $db->[0]->{$key} : undef;
}

sub STORE
{
	my $this = shift;
	my ($key, $value) = @_;
	my $new = 0;
	my $db = $this->fetch_db($key, \$new);
	$db->[0]->{$key} = $value;
	++$db->[1] if $new;
}

sub EXISTS
{
	my $this = shift;
	my ($key) = @_;
	return defined $this->find_db_by_key($key);
}

sub DELETE
{
	my $this = shift;
	my ($key) = @_;
	my $db = $this->find_db_by_key($key);
	if ($db) {
		delete $db->[0]->{$key};
		--$db->[1];
	}
}

sub CLEAR
{
	my $this = shift;
	foreach my $db (@{$this->{dbs}}) {
		foreach my $key (keys %{$db->[0]}) {
			delete $db->[0]->{$key};
		}
	}
}

sub FIRSTKEY
{
	my $this = shift;
	$this->{idx} = 0;
	return $this->NEXTKEY();
}

sub NEXTKEY
{
	my $this = shift;
	my $idx = \$this->{idx};
	my $db = $this->{dbs}[$$idx] || return undef;
	my ($key, $val) = each %{$db->[0]};
	if (!defined $key) {
		++$$idx;
		$db = $this->{dbs}[$$idx] || return undef;
		keys %{$db->[0]};
		($key, $val) = each %{$db->[0]};
	}
	return $key;
}

sub SCALAR
{
	my $this = shift;
	my $size = 0;
	foreach my $db (@{$this->{dbs}}) {
		$size += $db->[1];
	}
	return $size;
}

sub set_db_size
{
	my $this = shift;
	my ($new_db_size) = @_;
	die "Currently cannot resize databases when data has been loaded"
		if $this->SCALAR();
	$this->{db_size} = $new_db_size;
}

sub create_hash
{
	my $this = shift;
	my %hash = ();
	my $chunk = scalar(@{$this->{dbs}});
	my $filename = $this->{db_name} ?
		$this->{db_name} . ".$chunk" : undef;
	tie %hash, "DB_File", $filename, @{$this->{opts}};
	push @{$this->{dbs}}, [\%hash, 0];
}

sub init_hash
{
	my $this = shift;
	my $search_path = $this->{db_name} . ".*[0-9]";
	#try loading any existing databases
	while (my $existing_db = glob($search_path)) {
		my %hash;
		tie %hash, "DB_File", $existing_db, @{$this->{opts}};
		#store existing database plus its size
		push @{$this->{dbs}}, [\%hash, scalar(keys(%hash))];
	}
	#if no previous database found, create a new one
	$this->create_hash() if (!scalar(@{$this->{dbs}}));
}

sub fetch_db
{
	my $this = shift;
	my ($key, $new) = @_;
	my $db = $this->find_db_by_key($key);
	if ($db) {
		$$new = 0 if $new;
		return $db;
	}
	$$new = 1 if $new;
	$db = $this->{dbs}->[-1];
	$this->create_hash() if $db->[1] >= $this->{db_size};
	return $this->{dbs}->[-1];
}

sub find_db_by_key
{
	my $this = shift;
	my ($key) = @_;
	foreach my $db (@{$this->{dbs}}) {
		return $db if exists $db->[0]->{$key};
	}
	return undef;
}

1;
