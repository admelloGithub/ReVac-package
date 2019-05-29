package Coati::General;

use strict;
use File::Basename;

#################################

sub handle_btab_names {
    my ($self, $seq_id, $gene_id, $db) = @_;
	my $id;

    $self->{_logger}->debug("Args: ",join(',',@_)) if $self->{_logger}->is_debug;
    my ($btab_modif_date);
    my ($btab_file,$btab_dir, $id) = $self->{_backend}->get_handle_btab_names($seq_id, $gene_id, $db);
    if(!-e($btab_file) && -e($btab_file . ".gz")) { 
	$btab_modif_date = (stat($btab_file.".gz"))[9];
    } 
    else {
	$btab_modif_date = (stat($btab_file))[9];
    }
    my $btab_file_date = localtime($btab_modif_date);
    return ($btab_file, $btab_modif_date, $btab_file_date, $btab_dir, $id);
}

sub handle_blast_names {
    my ($self, $seq_id, $gene_id, $db) = @_;
    my ($blast_file, $blast_modif_date, $blast_file_date, $custom_blasts) = $self->{_backend}->get_handle_blast_names($seq_id, $gene_id, $db);
    return ($blast_file, $blast_modif_date, $blast_file_date, $custom_blasts);
}

sub handle_repeat_viewer_names {
    my ($self, $seq_id, $db) = @_;
    
    $self->{_logger}->debug("Args: ",join(',',@_)) if $self->{_logger}->is_debug;

    my $rv_file = sprintf("%s/%s/asmbls/%s/$seq_id.contig.repeatlist",
			  $ENV{ANNOTATION_DIR},
			  uc($db),
			  $seq_id,
			  $seq_id);
    return $rv_file;
}

sub run_signalp {
    my($self, $type, $trunc, $fasta_file) = @_;
    $self->{_logger}->debug("Args: ",join(',',@_)) if $self->{_logger}->is_debug;

    my $command = sprintf("$ENV{SIGNAL_P} -t '%s' -m 'nn+hmm' -trunc %s %s", $type, $trunc, $fasta_file);

    return $command;
}

sub revcomp_coord {
    my ($self, $coord, $seq_length) = @_;
    $self->{_logger}->debug("Args: ",join(',',@_)) if $self->{_logger}->is_debug;
    return ($seq_length - $coord + 1);
}

sub reverse_complement {
    my ($self, $seq_ref) = @_;
    $self->{_logger}->debug("Args: ",join(',',@_)) if $self->{_logger}->is_debug;
    my $rc = reverse($$seq_ref);
    $rc =~ tr/ACGTacgtyrkmYRKM/TGCAtgcarymkRYMK/;
    return \$rc;
}

sub get_orientation {
    my ($self, $end5, $end3) = @_;
    
    $self->{_logger}->debug("Args: ",join(',',@_)) if $self->{_logger}->is_debug;

    if ($end5 < $end3) {
        return '+';
    } 
    elsif ($end5 > $end3) {
        return '-';
    } 
    else {
        return undef;
    }
}

#################################

1;
