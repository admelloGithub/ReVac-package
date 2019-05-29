#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use warnings;
use Data::Dumper;
use XML::Twig;
use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;
use File::OpenFile qw(open_file);

my $in                  = "/dev/stdin";
my $out                 = *STDOUT;
my $output_mrna_feats   = 0;
my $extract_all_ec      = 0;
my $percent_n_cutoff    = 10;
my $feat_id_to_seq_id;
my $organism;
my $seq_id;
my %opts = ();
## genes with equal or more than this % of Ns will be skiped (warning printed)


parse_options();
convert();

sub print_usage
{
    my $progname = basename($0);
    die << "END";
usage: $progname [--input|-i <input_bsml>]
        [--output_dir|-o <output_tbl>] [--mrna|-m] [--ec_all|-e] [-h]

        -m: export mRNA features [ default - 0 (false) ]
        -e: extract all ec numbers (including incomplete e.g. #.#.#.-)
            [ default - false ]
END
}

sub parse_options
{
    GetOptions(\%opts, 
               "input|i=s", 
               "output_dir|o=s", 
               "mrna|m=i", 
               "help|h");
    print_usage() if $opts{help};
    $in = $opts{input} if $opts{input};

    my $file = basename($opts{input}, ('.bsml','.bsml.gz', '.bsml.gzip'));
    $out = open_file($opts{output_dir}."/$file.mugsymap", "out") or
        die "Error writing tbl to $opts{output}: $!"
        if $opts{output_dir};
    $output_mrna_feats = 0 unless $opts{mrna};
    $extract_all_ec = 1 if $opts{ec_all};
    if ( defined $opts{percent_n_cutoff} ) {
        $percent_n_cutoff = $opts{percent_n_cutoff};
    }


}

sub convert
{
    my $twig = new XML::Twig();
    my %feats = ();


    my $output_features = [];
    my $organism_name;

    # First pull out all of the features.
    $twig->setTwigRoots({
#        'Sequence[@class="assembly"]' => \&process_sequence,
        'Organism' => \&process_organism,
	'Sequence[@class="assembly"]' => \&process_sequence, 
        'Feature' => sub { process_feature(\%feats, @_); }});

    my $ifh = $in;
    if($opts{input}) {
        $ifh = open_file($in, 'in');
    }
    $twig->parse($ifh);
    $twig = new XML::Twig();
    close $ifh;
    # Next pull out all of the featuregroups
    $twig->setTwigRoots({'Feature-group' =>
            sub { process_feature_group(\%feats, $output_features, @_); } });
    if($opts{input}) {
        $ifh = open_file($in, 'in');
    }
    $twig->parse($ifh);
    close $ifh;

    # Sort the features based on location
    my @srted_feats = sort {$a->{'start'} <=> $b->{'start'}} @$output_features;
	if(@srted_feats == 0) {
			print $out "$seq_id|||-\t$seq_id\t-\t-\t-\t-\t-\t$organism\t-\n";
	}

    # Print them out 
    foreach my $feat (@srted_feats) {
        print $out join("\t", ("$feat->{'seq_id'}|||$feat->{'gene'}",$feat->{'seq_id'},$feat->{'start'},$feat->{'stop'},$feat->{'strand'},$feat->{'polypeptide_id'},$feat->{'gene_id'},$organism,$feat->{'gene_product'}))."\n";
    }
}

sub process_feature
{
    my ($feats, $twig, $elt) = @_;
    my $id = $elt->att('id');
    $feat_id_to_seq_id->{$id} = $elt->parent('Sequence')->att('id');
    $feats->{$id} = $elt;
    $twig->purge();
}

sub process_sequence
{
    my ($twig, $elt) = @_;
    $seq_id = $elt->att('id');
    print STDERR  "Have sequence $seq_id from $organism\n";
    $twig->purge();
}

sub process_organism
{
    my ($twig, $elt) = @_;
    my $id = $elt->att('id');
    my $genus = $elt->att('genus');
    my $spec = $elt->att('species');
    my $strain = $elt->first_child('Strain');
    if($strain) {
        $spec .= " ".($strain->first_child('Attribute')->att('content'));
    }

    my $org = $genus." ".$spec;
    $org =~ s/[\s\/]+/_/g;
    $org =~ s/[\.;]//g;
    $organism = $org;

    $twig->purge();
}

sub process_feature_group
{
    my ($feats, $output_features, $twig, $elt) = @_;
    my @exons = ();
    my $cds = undef;
    my $class = undef;
    my $gene = undef;
    my $polypeptide = undef;
    my $transcript = undef;
    my @repeats = ();
    my $group_by_class = {};

    foreach my $feat_member ($elt->children('Feature-group-member')) {
        my $feat_type = $feat_member->att('feature-type');
        my $featref = $feat_member->att('featref');
        my $feat = $feats->{$featref};
        if(!$group_by_class->{$feat_type}) {
            $group_by_class->{$feat_type} = [];
        }
        push(@{$group_by_class->{$feat_type}}, $feat);
        
        if ($feat_type eq "exon") {
            push @exons, $feat;
        }
        elsif ($feat_type eq "transcript") {
            $class = "transcript";
            $transcript = $feat;
        }
        elsif ($feat_type eq "tRNA") {
            $class = "tRNA";
            $transcript = $feat;
        }
        elsif ($feat_type eq "ncRNA") {
            $class = "ncRNA";
        }
        elsif ($feat_type eq "rRNA") {
            $class = "rRNA";
            $transcript = $feat;
        }
        elsif ($feat_type eq "gene") {
            $gene = $feat;
        }
        elsif ($feat_type eq "CDS") {
            $cds = $feat;
        }
        elsif ($feat_type eq "polypeptide") {
            $polypeptide = $feat;
        }
        elsif ($feat_type eq "repeat_region") {
            push @repeats, $feat;
        } elsif( $feat_type eq "signal_peptide" ) {
            $class = "signal_peptide";
        }
    }


    if(!$transcript) {
        print STDERR "Unable to find a transcript for group".$elt->att('id')."\n";
        return;
    }
    # Going to pull the location info here
    my $location = $transcript->first_child('Interval-loc');
    my $fmin = $location->att('startpos');
    my $fmax = $location->att('endpos');    
    my $strand = $location->att('complement') ==0 ? 1 : -1;
    my $aa_length = int(($fmax-$fmin)/3); # HACK - obviously this is not always right. The issue is that .ptt techincally wants amino acid length.

    # This is a HACK of course cause if there are more than one of these types we will have issues
    my $pid;

    foreach my $type (('transcript')) {
        my $pid_feat = $group_by_class->{$type}->[0];
        
        if($pid_feat) {
            # Need to be able to handle additional fields here.
            $pid = $pid_feat->att('id');
        }
        last if $pid;
    }
    
    # Pulling the gene product off of the transcript.
    my $att = $transcript->first_child('Attribute[@name="gene_product_name"]');
    my $gene_product;
    if(!$att) {
        $att = $gene->first_child('Attribute[@name="gene_product_name"]');
    }
    if($att) {
        $gene_product = $att->att('content');
    }
    else {
        print STDERR "Couldn't find gene product for $pid\n";
        $gene_product = 'None specified';
    }

    my $gene_val = '-';
    foreach my $type (('gene','transcript','CDS')) {
        my $gene_feat = $group_by_class->{$type}->[0];
        if($gene_feat) {
            my @xrefs = $gene_feat->children('Cross-reference');
            map {
                if(($_->att('identifier-type')) && ($_->att('identifier-type') eq 'locus') && $_->att('identifier') ne "") {
                    $gene_val = $_->att('identifier');
                }
                elsif(($_->att('database')) && ($_->att('database') eq 'NCBILocus')) {
                    $gene_val = $_->att('identifier');
                }
                elsif(($_->att('database')) && ($_->att('database') eq 'NCBI_locus_tag')) {
                    $gene_val = $_->att('identifier');
                }
                elsif(($_->att('database')) && ($_->att('database') eq 'NCBI_locus_tag_public_locus')) {
                    $gene_val = $_->att('identifier');
                }
                elsif(($_->att('database')) && ($_->att('database') eq 'NCBI_gi')) {
                    $gene_val = $_->att('identifier');
                }
            }@xrefs;
        }
        last if $gene_val ne '-';
    }
    die "Unable to set name for gene: ",$gene->att('id') if($gene_val eq '-' || $gene_val eq '');
    if($pid) {
        push(@$output_features, {'start'          => $fmin,
                                 'stop'           => $fmax,
                                 'strand'         => $strand,
                                 'length'         => $aa_length,
                                 'pid'            => $pid,
                                 'gene'           => $gene_val,
                                 'gene_product'   => $gene_product,
                                 'polypeptide_id' => $polypeptide->att('id'),
                                 'gene_id'        => $gene->att('id'),
                                 'seq_id'         => $feat_id_to_seq_id->{$gene->att('id')}
             });
    }
}

sub process_misc_feats
{
    my ($feats) = @_;
    while (my ($key, $val) = each %{$feats}) {
        my $class = $val->att('class');
        if ($class eq "repeat_region") {
            print_feat([$val], "repeat_region");
        }
    }
}
