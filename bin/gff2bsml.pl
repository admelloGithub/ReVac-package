#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

# NOTES
#  -allow multiple transcripts/gene (each gets its own BSML feature group)
#  -allow id-less CDS features (with --use_cds_parent_ids), so long as parent ids are unique
#  -removed analysis and analysis links
#  -added miRNA support (and, more generally, support for multiple transcript-level types)
#  -ditto for tRNA, ncRNA, snoRNA, snRNA, rRNA
#  -made $feat_type_map into a global variable
#  -added --default_seq_class
#  -added --peptide_fasta, --peptide_id_regex, --peptide_id_prefix, --peptide_id_suffix to handle external file of peptide seqs
#  -added --dna_fasta, --dna_id_regex, --dna_id_prefix, --dna_id_suffix to handle external file of DNA seqs
#  -removed "temp_id" strings from feature element titles
#  -added 'organism' input that gets converted into a BSML <Genome> section
#  -added --organism_genetic_code, --organism_mt_genetic_code
#   -all Sequences link to that <Genome>
#   -reference Sequence lengths are determined from the GFF, if possible
#  -modified Seq-data-imports to use absolute paths
#  -added --insert_missing_mrnas, --insert_missing_exons
#  -added --feat_xref_mappings, --seq_xref_mappings

# Current restrictions/assumptions
#  -all sequences/features belong to 1 genome, the one named by --organism param

# TODO
#  -generalize to regular GFF?
#  -add a strict mode in which deviations from the GFF3/canonical gene spec. are reported
#   -see http://www.sequenceontology.org/gff3.shtml
#  -if --use_cds_parent_ids is in effect the script should check that the CDS features under each mRNA are
#   nonoverlapping: if not it may be impossible to correctly reconstruct the original gene model without
#   correctly-assigned CDS ids
#  -add options to control the handling of '*' characters in the peptide sequences
#   -current legacy2bsml/genbank2bsml scripts appear to be inconsistent wrt whether '*' is included
#   -and if the '*' is included then it gets counted in the seqlen in chado
#  -create Research/Analysis elements for computed elements in the BSML (add gene_model analysis param?)
#  -add --ref_seq_type param? change name of --default_seq_class param
#  -need a way to specify which features get linked to Analyses in the <Research> section
#  -replace die with logdie where permissible
#  -make --insert_polypeptides work in the presence of alternate splicing
#   (so long as --inserted_polypeptide_id_type=CDS)

=head1 NAME

gff32bsml.pl - Convert GFF3 annotation data into BSML.

=head1 SYNOPSIS

gff32bsml.pl
         --input=/path/to/annotation.gff3
         --output=/path/to/results.bsml
         --project=rca1
         --id_repository=/usr/local/projects/rca1/workflow/project_id_repository
         --organism='Aspergillus nidulans'
        [--peptide_fasta=/path/to/peptide-seqs.fsa
         --peptide_id_regex='^>(\S+)'
         --peptide_id_prefix=''
         --peptide_id_suffix='-Protein'
         --peptide_no_seqdata
         --dna_fasta=/path/to/dna-seqs.fsa
         --dna_id_regex='^>(\S+)'
         --dna_id_prefix='Chr'
         --dna_id_suffix=''
         --dna_no_seqdata
         --gene_type=gene
         --organism_genetic_code=1
         --organism_mt_genetic_code=4
         --use_cds_parent_ids
         --insert_missing_mrnas
         --insert_missing_exons
         --insert_polypeptides
         --inserted_polypeptide_id_type=CDS
         --default_seq_class=assembly
         --default_seq_attributes='SO:contig'
         --feat_attribute_mappings='Name:gene_product_name,comment:comment'
         --clone_feat_attribute_mappings='gene:gene_product_name:transcript,gene:gene_product_name:CDS'
         --feat_xref_mappings='ID:mydb:accession,Name:mydb:pub_locus'
         --clone_feat_xref_mappings='gene:NCBILocus:transcript,gene:NCBILocus:CDS'
         --seq_xref_mappings='ID:mydb:seq_id'
         --atts_with_unescaped_commas='ID,Name'
         --allow_genes_with_no_cds
         --sofa_type_mappings='transcript=>CDS,ORF=>gene'
         --log=/path/to/some.log
         --debug=4
         --help
         --man ]

=head1 OPTIONS

B<--input, -g>
    path to the input GFF file.

B<--peptide_fasta, -e>
    optional.  path to a multi-FASTA file containing the polypeptides for the genes in the GFF file.
    if this argument is supplied then the named file _must_ contain all the sequences for the 
    polypeptide/protein features in the input GFF file, _except_ those that appear in the ##FASTA
    section of the GFF file.

B<--peptide_id_regex>
    optional.  regular expression used to parse peptide unique id from FASTA deflines in --peptide_fasta
    and the FASTA section of the GFF file itself.  in the latter case, the regular expression must NOT
    match any non-peptide ids, otherwise --peptide_id_prefix and --peptide_id_suffix may be misapplied.
    default = ^>(\S+)

B<--peptide_id_prefix>
    optional. prefix to prepend to the id parsed by --peptide_id_regex to produce the polypeptide ID 
    used by the corresponding GFF feature.

B<--peptide_id_suffix>
    optional. suffix to append to the id parsed by --peptide_id_regex to produce the polypeptide ID 
    used by the corresponding GFF feature.

B<--peptide_no_seqdata>
    do not create Seq-data or Seq-data-import entries for the polypeptide sequences; use --peptide_fasta,
    if specified, only to determine the sequence length(s).

B<--dna_fasta>
    optional.  path to a multi-FASTA file containing the DNA sequences for the reference sequences
    in the GFF file.  if this argument is supplied then the named file _must_ contain all of the 
    sequences for the genomic sequence features in the input GFF file, _except_ those that appear
    in the ##FASTA section of the GFF file.

B<--dna_id_regex>
    optional.  regular expression used to parse the sequence unique id from FASTA deflines in 
    --dna_fasta and the FASTA section of the GFF file itself.  in the latter case, the regular 
    expression must NOT match any non-peptide ids, otherwise --dna_id_prefix and --dna_id_suffix 
    may be misapplied..   default = ^>(\S+)

B<--dna_id_prefix>
    optional. prefix to prepend to the id parsed by --dna_id_regex to produce the sequence ID 
    used by the corresponding GFF feature.

B<--dna_id_suffix>
    optional. suffix to append to the id parsed by --dna_id_regex to produce the sequence ID 
    used by the corresponding GFF feature.

B<--dna_no_seqdata>
    do not create Seq-data or Seq-data-import entries for the DNA/reference sequences; use --dna_fasta,
    if specified, or the sequences present in the GFF file, only to determine the sequence length(s).

B<--gene_type>
    optional.  specifies a SO type in the input GFF file that should be replaced with "gene".  this can
    be used to handle nonstandard inputs in which, for example, the term "match" is used as the parent
    of one or more CDS features (e.g., as in GeneWise GFF output)  Note that the SO types _are_
    case-sensitive.

B<--organism_genetic_code>
    optional.  specifies a genetic code value to appear in the BSML <Organism> element for --organism
    see http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c#SG1 for values

B<--organism_mt_genetic_code>
    optional.  specifies a mitochondrial genetic code value to appear in the BSML <Organism> element 
    for --organism.  see http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c#SG1 for values

B<--output, -o>
    output BSML file (which should not exist before running the script.)

B<--project, -p>
    project/database name used to create unique feature identifiers.

B<--id_repository, -i>
    path to the Ergatis project_id_repository for the current project.

B<--organism, -r>
    organism to place in the BSML <Genomes> section.  the first word will be used as the BSML genus,
    and the rest of the string will be placed in the BSML species attribute.

B<--use_cds_parent_ids,-c>
    optional.  whether to use the GFF parent id to identify CDS features that lack an ID of their own.

B<--insert_missing_mrnas>
    optional.  automatically insert missing mRNA features wherever they are missing from the GFF 
    (i.e., anytime a CDS feature is found hanging directly off a gene feature)
    
B<--insert_missing_exons>
    optional.  automatically inserts missing exon features in any gene model that does not have _any_
    exons but does have CDS features.  in this case a single exon will be created to exactly span each
    CDS.

B<--insert_polypeptides>
    optional.  automatically inserts a polypeptide feature in any gene model that has one or more CDS
    features but no polypeptide feature.

B<--inserted_polypeptide_id_type>
    optional. specifies which feature ('CDS', 'gene', 'transcript', or 'mRNA') the newly-inserted
    polypeptide ids should be based on (default = 'gene')  the id of the new polypeptide will be 
    formed by taking the GFF id of the corresponding named feature and appending "-Protein"  
    Note that "transcript" and "mRNA" are synonymous for the purposes of this option.

B<--default_seq_class,-e>
    optional.  default sequence class/SO type (e.g., 'assembly', 'supercontig') to use for sequences 
    for which the type cannot be parsed from the GFF file.

B<--default_seq_attributes,-e>
    optional.  comma-delimited list of BSML attributes (in colon-delimited key:value form) to associate
    with each genomic sequence in the GFF file.

B<--feat_attribute_mappings>
    optional.  a comma-delimited list of GFF attribute -> BSML Attribute mappings, each of which is
    defined by a source (GFF) attribute name and a target (BSML) attribute name, separated by a colon,
    for example "Name:gene_product_name" will take each BSML 'Name' value and insert it into the 
    BSML document as an <Attribute> with name="Name" and content=the corresponding GFF attribute value.
    Any embedded colons may be escaped with a backslash ("\:")

B<--clone_feat_attribute_mappings>
    optional.  a comma-delimited list of colon-separated values.  each member of the list is a colon-
    separated triplet that contains 1. a source SOFA type, 2. a BSML attribute name, and 3. a target 
    SOFA type.  e.g., "gene:gene_product_name:CDS,gene:gene_product_name:transcript"  This particular
    example specifies that--within each gene model/feature group--any gene_product_name BSML Attributes
    attached to the gene feature will be copied over to the associated transcript and CDS features
    (without introducing duplicate attribute values).  This cloning/copying is done _after_ any applicable 
    --feat_attribute_mappings have been processed for all of the involved features.

B<--feat_xref_mappings,-x>
    optional.  a comma-delimited list of GFF attribute -> BSML Cross-reference mappings, each of which
    is defined by a set of 3 values separated by colons, as in 'ID:mydb:accession', which specifies
    that each BSML 'ID' value should be inserted into the BSML document as a <Cross-reference> 
    element with database="mydb" and identifier-type="accession" (and identifier set to the ID value.)
    Any embedded colons may be escaped with a backslash ("\:")
  
B<--clone_feat_xref_mappings>
    optional.  behaves similarly to --clone_feat_attribute_mappings, but copying BSML <Cross-references>
    instead of <Attributes>

B<--seq_xref_mappings>
    optional.  identical to --xref_mappings except that it defines cross-reference mappings for the 
    reference _sequences_ in the GFF files (i.e., those things whose IDs appear in GFF column 1)
    in the case where the sequences are not themselves defined as GFF features the only attribute for
    which a mapping may be specified is 'ID'.

B<--atts_with_unescaped_commas>
    optional.  a comma-delimited list of GFF3 attributes whose values (incorrectly) contain unescaped
    commas.  instead of parsing these values as comma-delimited lists the program will treat each as
    a single value and will parse it as though the commas embedded in the value had been correctly 
    URL-escaped.

B<--allow_genes_with_no_cds>
    optional.  normally the script expects that any non-protein-coding gene will either: 1. be assigned
    an explicit SO type of 'pseudogene' instead of 'gene' or 2. have a non-coding RNA SO type at the
    transcript level (e.g., ncRNA, tRNA, etc.)  therefore if the script is running with --insert_polypeptides
    and finds a gene for which neither of these things are true, it will fail with an error if the gene
    lacks a CDS.  use this flag to permit the conversion to continue with a warning instead.

B<--sofa_type_mappings>
    optional.  a comma-delimited list of SOFA id or name mappings, in the form id1=>id2 or name1=>name2.
    each occurrence of id1 or name1 in the GFF3 SOFA feature type column (GFF3 column #3) will be
    replaced by id2/name2 prior to performing any subsequent parsing.  this option can be useful for
    making a GFF3 file more closely match the canonical gene encoding.  for example, if a GFF3 file
    uses the term "ORF" instead of "gene" it can be rectified by specifying
    --sofa_type_mappings='ORF=>gene'

B<--log,-l>
    optional.  path to a log file the script should create.  will be overwritten if
    it already exists.

B<--debug,-d>
    optional.  the debug level for the logger (an integer)

B<--help,-h>
    print usage/help documentation

B<--man,-m>
    print detailed usage/help documentation

=head1 DESCRIPTION

This script parses a single GFF input file containing annotation data (for one or more
genomic sequences) and writes it out as a BSML document that can be subsequently fed 
into bsml2chado to load it into a Chado database.  This script was originally created 
by copying and then generalizing Brett Whitty's evmgff32bsml.pl utility.

=head1 INPUT

A single GFF input file.  See the GFF documentation for more details.  Note that the
script has a number of options that allow it to correctly parse/convert a variety of
GFF-format files, even those that may deviate from the official specification in one
way or another.

=head1 OUTPUT

A single BSML document.  See the BSML documentation for more details.

=head1 CONTACT

    Jonathan Crabtree
    jonathancrabtree@gmail.com

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Data::Dumper;
use BSML::BsmlBuilder;
use Ergatis::IdGenerator;
use Ergatis::Logger;
use FileHandle;
use File::Spec;
use URI::Escape;

## global/default values    
my $DEFAULT_SEQ_CLASS = 'assembly';
my $DEFAULT_PEPTIDE_ID_REGEX = '^>(\S+)';
my $DEFAULT_DNA_ID_REGEX = '^>(\S+)';
my $GFF_SEQID_REGEX = '^(>\S.*)';
my $GLOBAL_ID_COUNTER = 0;
my $NODES = {};
# mapping from GFF sequence id to BSML sequence id.  these need not be the same.
my $GFF2BSML_SEQID_MAP = {};

# mapping from GFF3 feature type to Chado cvterm.  if the script encounters
# a GFF type not found in this hashref it will fail
my $FEAT_TYPE_MAP = {
    'gene' => 'gene',
    'ORF' => 'gene',
    'CDS' => 'CDS',
    'cds' => 'CDS',
    'exon' => 'exon',
    'intron' => 'intron',
    'mRNA' => 'transcript',               
    'tmRNA' => 'tmRNA',               
    'transcript' => 'transcript',               
    'miRNA' => 'miRNA',
    'protein' => 'polypeptide',
    'polypeptide' => 'polypeptide',
    'tRNA' => 'tRNA',
    'ncRNA' => 'ncRNA',
    'snoRNA' => 'snoRNA',
    'snRNA' => 'snRNA',
    'rRNA' => 'rRNA',
    'three_prime_UTR' => 'three_prime_UTR',
    'five_prime_UTR' => 'five_prime_UTR',
    'three_prime_utr' => 'three_prime_UTR',
    'five_prime_utr' => 'five_prime_UTR',
    'gap' => 'gap',
    'pseudogene' => 'pseudogene',
    'contig' => 'contig',
    'supercontig' => 'supercontig',
    'chromosome' => 'chromosome',
};

# list of types that are allowed for the mRNA; must be a subset of those
# listed _as values_ in $FEAT_TYPE_MAP
my $MRNA_TYPES = {
    'transcript' => 1,
    'miRNA' => 1,
    'tRNA' => 1,
    'ncRNA' => 1,
    'snoRNA' => 1,
    'snRNA' => 1,
    'rRNA' => 1,
    'tmRNA' => 1,
};
# TODO - add all posible mRNA types from SO?

my $NONCODING_RNA_TYPES = {
    'miRNA' => 1,
    'tRNA' => 1,
    'ncRNA' => 1,
    'snoRNA' => 1,
    'snRNA' => 1,
    'rRNA' => 1,
    # not exactly true, but it wouldn't be handled correctly in the other category:
    'tmRNA' => 1,
};

## input/options
my $options = {};
my $results = GetOptions($options, 
                         'input|g=s',
                         'peptide_fasta=s',
                         'peptide_id_regex=s',
                         'peptide_id_prefix=s',
                         'peptide_id_suffix=s',
                         'peptide_no_seqdata',
                         'dna_fasta=s',
                         'dna_id_regex=s',
                         'dna_id_prefix=s',
                         'dna_id_suffix=s',
                         'dna_no_seqdata',
                         'gene_type=s',
                         'organism_genetic_code=i',
                         'organism_mt_genetic_code=i',
                         'output|o=s',
                         'project|p=s',
                         'id_repository|i=s',
                         'organism|r=s',
                         'use_cds_parent_ids|c',
                         'insert_missing_mrnas',
                         'insert_missing_exons',
                         'insert_polypeptides',
                         'inserted_polypeptide_id_type=s',
                         'default_seq_class|e=s',
                         'default_seq_attributes|e=s',
                         'feat_attribute_mappings=s',
                         'clone_feat_attribute_mappings=s',
                         'feat_xref_mappings|x=s',
                         'clone_feat_xref_mappings|x=s',
                         'seq_xref_mappings=s',
                         'atts_with_unescaped_commas=s',
                         'allow_genes_with_no_cds',
                         'sofa_type_mappings=s',
                         'log|l=s',
                         'debug|d=i',
                         'help|h',
                         'man|m'
                         ) || pod2usage();

&pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT}) if ($options->{'help'});
&pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if ($options->{'man'});

&check_parameters($options);

my $unescaped_comma_atts = {};
if (defined($options->{'atts_with_unescaped_commas'})) {
    map { $unescaped_comma_atts->{$_} = 1; } split(/\s*,\s*/, $options->{'atts_with_unescaped_commas'});
}

if (defined($options->{'gene_type'})) {
    $FEAT_TYPE_MAP->{$options->{'gene_type'}} = $FEAT_TYPE_MAP->{'gene'};
}

## initialize logging
my $logfile = $options->{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile, 'LOG_LEVEL'=>$options->{'debug'});
$logger = Ergatis::Logger::get_logger();

if (defined($options->{'sofa_type_mappings'})) {
    my @mappings = split(/,/, $options->{'sofa_type_mappings'});
    foreach my $mapping (@mappings) {
        my($from,$to) = ($mapping =~ /^(\S+)\s*\=\>\s*(\S+)$/);
        $logger->logdie("unable to parse sofa_type_mapping from $mapping") if (!defined($to));
        $FEAT_TYPE_MAP->{$from} = $to;
        $logger->debug("sofa_type_mapping: $from => $to");
    }
}

## MAIN SECTION
my $idcreator = Ergatis::IdGenerator->new('id_repository' => $options->{'id_repository'});
my @root_nodes = ();
my($genus, $species) = ($options->{'organism'} =~ /^(\S+) (.*)$/);

## parse peptide FASTA (optional)
my $peptides = undef;
my $pep_prefix = $options->{'peptide_id_prefix'};
my $pep_suffix = $options->{'peptide_id_suffix'};

if (defined($options->{'peptide_fasta'})) {
    my $pep_id_fn = sub {
        my($id) = @_;
        my $new_id = (defined($pep_prefix)) ? $pep_prefix . $id : $id;
        $new_id .= $pep_suffix if (defined($pep_suffix));
        return $new_id;
    };
    $peptides = &read_sequence_lengths($options->{'peptide_fasta'}, $options->{'peptide_id_regex'}, $pep_id_fn);
}

## parse genomic DNA FASTA (optional)
my $dna_seqs = undef;
my $dna_prefix = $options->{'dna_id_prefix'};
my $dna_suffix = $options->{'dna_id_suffix'};

if (defined($options->{'dna_fasta'})) {
    my $dna_id_fn = sub {
        my($id) = @_;
        my $new_id = (defined($dna_prefix)) ? $dna_prefix . $id : $id;
        $new_id .= $dna_suffix if (defined($dna_suffix));
        return $new_id;
    };
    $dna_seqs = &read_sequence_lengths($options->{'dna_fasta'}, $options->{'dna_id_regex'}, $dna_id_fn);
}

## parse GFF
open (IN, $options->{'input'}) || die "couldn't open input file '$options->{input}'";
my $lnum = 0;
my $found_fasta_header = 0;

while (<IN>) {
    chomp;
    ++$lnum;
    
    # fasta-format sequences
    if (/^>/ || /^\#.*fasta.*/i) {
        $found_fasta_header = 1;
        last;
    }
    # ignore blank lines and those beginning with '#'
    elsif (/^\#|^\s*$/) {
        # no-op
    } 
    else {
        my $record = parse_record($_, $lnum);

        # if it has no parents it's a root node    
        if ((!$record->{'Parent'}) && (!$record->{'Derives_from'})) {
            push(@root_nodes, $record);
        }
    
        # store records in a set of arrays indexed by ID
        if ($record->{'ID'}) {
            push(@{$NODES->{$record->{'ID'}}->{'records'}},$record);
        }
    }
}

# write FASTA sequences from GFF directly to an external file
my $gff_seqs = undef;
my $gff_seq_file = undef;

if ($found_fasta_header) {
    $gff_seq_file = $options->{'output'} . ".fsa";
    $logger->logdie("$gff_seq_file already exists; please remove it and rerun") if (-e $gff_seq_file);
    my $ffh = FileHandle->new();
    $ffh->open(">$gff_seq_file") || die "";
    while (my $line = <IN>) {
        $ffh->print($line);
    }
    $ffh->close();

    my $pep_regex = $options->{'peptide_id_regex'};
    my $dna_regex = $options->{'dna_id_regex'};

    # then parse that file back to get the sequences
    my $gff_seqid_fn = sub {
        my($id) = @_;

        if (($pep_regex ne $DEFAULT_PEPTIDE_ID_REGEX) && ($id =~ /$pep_regex/)) {
            my $new_id = (defined($pep_prefix)) ? $pep_prefix . $1 : $1;
            $new_id .= $pep_suffix if (defined($pep_suffix));
            return $new_id;    
        }
        elsif (($dna_regex ne $DEFAULT_DNA_ID_REGEX) && ($id =~ /$dna_regex/)) {
            my $new_id = (defined($dna_prefix)) ? $dna_prefix . $1 : $1;
            $new_id .= $dna_suffix if (defined($dna_suffix));
            return $new_id;    
        }

        # default
        my ($new_id) = ($id =~ /^>(\S+)/);
        return $new_id; 
    };
    $gff_seqs = &read_sequence_lengths($gff_seq_file, $GFF_SEQID_REGEX, $gff_seqid_fn);
}

close IN;

# populate children for each record
foreach my $id(keys(%{$NODES})) {
    foreach my $record(@{$NODES->{$id}->{'records'}}) {
        if (defined($record->{'children'})) {
            $logger->warn("about to overwrite children field of record $record with id=$id)");
        }
        $record->{'children'} = $NODES->{$id}->{'children'};
        delete $NODES->{$id}->{'children'};
    }
}

my $gene_nodes = [];
fetch_node_type('gene', \@root_nodes, $gene_nodes);
my $num_genes = scalar(@$gene_nodes);
$logger->debug("found $num_genes gene(s)");
my $doc = new BSML::BsmlBuilder();
my $seq2feat_table = undef; # global var

# add Genomes
my $genome = $doc->createAndAddGenome();
my $genome_id = $genome->{attr}->{id};
my $organism = $doc->createAndAddOrganism('genome' => $genome, 'genus' => $genus, 'species' => $species);
my $code = $options->{'organism_genetic_code'};
my $mt_code = $options->{'organism_mt_genetic_code'};

my $code_att = $doc->createAndAddBsmlAttribute($organism, 'genetic_code', $code) if (defined($code));
my $mt_code_att = $doc->createAndAddBsmlAttribute($organism, 'mt_genetic_code', $mt_code) if (defined($mt_code));

# copied this from legacy2bsml.pl.  not clear why this isn't handled internally to the BsmlDoc 
# (or why 'xrefctr' is used for an id shared by both x-refs and genomes)
$doc->{'xrefctr'}++;

# database cross-reference mappings, indexed by GFF attribute
my $feat_xref_mappings = &parse_xref_mappings($options->{'feat_xref_mappings'});
my $seq_xref_mappings = &parse_xref_mappings($options->{'seq_xref_mappings'});

# GFF -> BSML Attribute mappings
my $feat_att_mappings = &parse_attribute_mappings($options->{'feat_attribute_mappings'});

# feature xrefs and attributes that must be cloned
my $clone_feat_att_mappings = &parse_clone_attribute_mappings($options->{'clone_feat_attribute_mappings'});
my $clone_feat_xref_mappings = &parse_clone_xref_mappings($options->{'clone_feat_xref_mappings'});

foreach my $root (@{$gene_nodes}) {
    if ($root->{'_type'} eq 'gene') {

        # insert missing features if the appropriate options are set
        my ($cds_type, $mrna_type, $exon_type, $intron_type, $utr3_type, $utr5_type) = map {$FEAT_TYPE_MAP->{$_}} ('CDS', 'mRNA', 'exon', 'intron', 'three_prime_UTR', 'five_prime_UTR');
        
        # insert mRNA, if none and --insert_missing_mrnas set
        if ((!defined($root->{'children'}->{$mrna_type}) && $options->{'insert_missing_mrnas'})) {
            $logger->debug("inserting missing mRNA for " . $root->{'ID'});
            my $trans_id = $root->{'ID'} . "-mRNA";
            my $cds_list = $root->{'children'}->{$cds_type};
            my $exon_list = $root->{'children'}->{$exon_type};
            my $intron_list = $root->{'children'}->{$intron_type};
            my $mrna_list = $root->{'children'}->{$mrna_type} = [];
            my $utr3_list = $root->{'children'}->{$utr3_type};
            my $utr5_list = $root->{'children'}->{$utr5_type};

            # extent of mRNA should = extent of exons (if present), of CDS feats if not
            my $feat_list = defined($exon_list) ? $exon_list : $cds_list;
            my($m_start, $m_end, $m_strand) = &union_feat_intervals($feat_list);

            if (!defined($exon_list) && !defined($cds_list)) {
                $logger->warn("Trying to insert mRNA $trans_id for gene " . $root->{'ID'} . " but no CDS or exons could be found: using gene coords instead.");
                ($m_start, $m_end, $m_strand) = ($root->{_start}, $root->{_end}, $root->{_strand});
            } elsif (($m_start != $root->{_start}) || ($m_end != $root->{_end}) || ($m_strand != $root->{_strand})) {
                $logger->debug($root->{'ID'} . ": mrna coords=$m_start-$m_end/$m_strand gene coords=" . $root->{_start} . "-" . $root->{_end} . "/" . $root->{_strand} . "\n");
            }
            
            push(@$mrna_list, {
                'ID' => $trans_id,
                '_seqid' => $root->{'_seqid'},
                '_type' => $mrna_type,
                'children' => {},
                'Parent' => [ $root->{'ID'} ],
                '_source' => $root->{'_source'},
                '_score' => '.',
                '_strand' => $m_strand,
                '_start' => $m_start,
                '_end' => $m_end,
                # TODO - propagation of cross-references to automatically-inserted features should be optional:
                '_record' => $feat_list->[0]->{'_record'},
            });
            
            # make mRNA the parent of the CDS, exon, and intron features
            map { push(@{$_->{'Parent'}}, $trans_id); } @$cds_list if (defined($cds_list));
            map { push(@{$_->{'Parent'}}, $trans_id); } @$exon_list if (defined($exon_list));
            map { push(@{$_->{'Parent'}}, $trans_id); } @$intron_list if (defined($intron_list));
            map { push(@{$_->{'Parent'}}, $trans_id); } @$utr3_list if (defined($utr3_list));
            map { push(@{$_->{'Parent'}}, $trans_id); } @$utr5_list if (defined($utr5_list));
        }

        # insert exons, if none and --insert_missing_exons set
        my $exon_list = $root->{'children'}->{$exon_type};
        my $num_exons = defined($exon_list) ? scalar(@$exon_list) : 0;
        if (($num_exons == 0) && $options->{'insert_missing_exons'}) {
            my $cds_list = $root->{'children'}->{$cds_type};
            my $mrna_list = $root->{'children'}->{$mrna_type};
            my $nm = scalar(@$mrna_list);

            if (defined($cds_list)) {
                if (scalar(@$mrna_list) != 1) {
                    $logger->error("can't insert exons for " . $root->{'ID'} . ": gene has $nm mRNAs");
                } else {
                    my $exon_list = $root->{'children'}->{$exon_type} = [];
                    my $nc = scalar(@$cds_list);
                    my $mrna_id = $mrna_list->[0]->{'ID'};
                    $logger->debug("inserting $nc missing exon(s) for " . $root->{'ID'});
                    my $mrna_exons = $mrna_list->[0]->{'children'}->{$exon_type};
                    if (!defined($mrna_exons)) {
                        $mrna_exons = $mrna_list->[0]->{'children'}->{$exon_type} = [];
                    }

                    my $enum = 1;
                    foreach my $cds (@$cds_list) {
                        my $exon_id = $mrna_id . "-e" . $enum++;
                        my $exon = {
                            'ID' => $exon_id,
                            '_seqid' => $cds->{'_seqid'},
                            '_type' => $exon_type,
                            'children' => {},
                            'Parent' => [ $mrna_id ],
                            '_source' => $cds->{'_source'},
                            '_score' => '.',
                            '_strand' => $cds->{'_strand'},
                            '_start' => $cds->{'_start'},
                            '_end' => $cds->{'_end'},
                            # TODO - propagation of cross-references to automatically-inserted features should be optional:
                            '_record' => $cds->{'_record'},
                        };
                        push(@$exon_list, $exon);
                        push(@$mrna_exons, $exon);
                    }
                }
            } else {
                $logger->debug("not inserting exons for $root->{ID}: no CDS seqs found");
            }
        }

        my $features = {};
        process_node($root, $features);
        gene_feature_hash_to_bsml($NODES, $features, $gff_seqs, $peptides, $dna_seqs, $root->{'_seqid'});
    }
}

# add the analysis element
# TODO - need a way to specify which features get linked to an analysis
#$doc->createAndAddAnalysis(
#                            id => 'EVM_analysis',
#                            sourcename => $options->{'output'},
#                          );

$doc->write($options->{'output'});
exit(0);

## subroutines

sub parse_xref_mappings {
    my($mapping_str) = @_;
    return &parse_mappings($mapping_str, 'xref', ['db_name', 'id_type']);
}

sub parse_attribute_mappings {
    my($mapping_str) = @_;
    return &parse_mappings($mapping_str, 'attribute', ['bsml_att']);
}

sub parse_clone_attribute_mappings {
    my($mapping_str) = @_;
    return &parse_mappings($mapping_str, 'cloned attribute', ['bsml_att', 'target_type']);
}

sub parse_clone_xref_mappings {
    my($mapping_str) = @_;
    return &parse_mappings($mapping_str, 'cloned xrefs', ['db_name', 'target_type']);
}

sub parse_mappings {
    my($mapping_str, $type, $fieldnames) = @_;
    my $mappings = {};
    my $nfn = scalar(@$fieldnames);

    if (defined($mapping_str)) {
        my @mstrs = split(',', $mapping_str);
        foreach my $mstr (@mstrs) {
            my($key, @fields) = split(':', $mstr);
            $key =~ s/\://g;
            my $list = $mappings->{$key};
            $list = $mappings->{$key} = [] if (!defined($list));
            my $val = {};
            my $ff = scalar(@fields);
            die "expected $nfn fields, only found $ff in '$mstr'" if ($nfn != scalar(@fields));

            for (my $i = 0;$i < $nfn;++$i) {
                $val->{$fieldnames->[$i]} = $fields[$i];
            }
            push(@$list, $val);
            $logger->info("parsed $type mapping for key: $key, " . join(", ", map {$_ . "=" . $val->{$_}} @$fieldnames));
        }
    }
    return $mappings;
}

sub union_feat_intervals {
    my($feats) = @_;
    my $min_start = undef;
    my $max_end = undef;
    my $strand = undef;
    
    foreach my $feat (@$feats) {
        my($fs, $fe, $fstr) = map { $feat->{$_} } ('_start', '_end', '_strand');
        die "feature $feat->{ID} has no _start" if (!defined($fs));
        die "feature $feat->{ID} has no _end" if (!defined($fe));
        die "feature end < start in feature $feat->{ID}" if ($fe < $fs);
        die "feature $feat->{ID} has no _strand" if (!defined($fstr));

        if (!defined($min_start)) {
            $min_start = $fs;
            $max_end = $fe;
            $strand = $fstr;
        } else {
            die "strand mismatch ($strand vs $fstr) in feature $feat->{ID}" if ($strand != $fstr);
            $min_start = $fs if ($fs < $min_start);
            $max_end = $fe if ($fe > $max_end);
        }
    }

    return($min_start, $max_end, $strand);
}

sub check_parameters {
	my $options = shift;

    ## make sure required parameters were passed
    my @required = qw(input output project id_repository organism);
    for my $option ( @required ) {
        unless  ( defined $options->{$option} ) {
            die "--$option is a required option";
        }
    }

    ## default values
    $options->{'default_seq_class'} = $DEFAULT_SEQ_CLASS if (!defined($options->{'default_seq_class'}));
    $options->{'project'} =~ tr/A-Z/a-z/;
    $options->{'dna_id_regex'} = $DEFAULT_DNA_ID_REGEX if (!defined($options->{'dna_id_regex'}));
    $options->{'peptide_id_regex'} = $DEFAULT_PEPTIDE_ID_REGEX if (!defined($options->{'peptide_id_regex'}));
    $options->{'inserted_polypeptide_id_type'} = 'gene' if (!defined($options->{'inserted_polypeptide_id_type'}));

    ## check files
    my $input = $options->{'input'};
    die "input GFF file ($input) does not exist" if (!-e $input);
    die "input GFF file ($input) is not readable" if (!-r $input);

    my $peps = $options->{'peptide_fasta'};
    if (defined($peps)) {
        die "input peptide FASTA file ($peps) does not exist" if (!-e $peps);
        die "input peptide FASTA file ($peps) is not readable" if (!-r $peps);
    }

    my $dnas = $options->{'dna_fasta'};
    if (defined($dnas)) {
        die "input DNA sequence FASTA file ($dnas) does not exist" if (!-e $dnas);
        die "input DNA sequence FASTA file ($dnas) is not readable" if (!-r $dnas);
    }

    if ($options->{'organism'} !~ /^\S+\s\S+/) {
        die "--organism must specify a species and genus, e.g., --organism='Aspergillus nidulans'";
    }

    if ($options->{'inserted_polypeptide_id_type'} !~ /^gene|CDS|mRNA|transcript$/) {
        die "--inserted_polypeptide_id_type must be one of the following: gene,CDS,mRNA,transcript";
    }
}

# Parse a single GFF3 record.
#
# $line - A chomp'ed GFF3 feature line.
#
sub parse_record {
    my ($line, $linenum) = @_;

    my $record = {};
    my $attrib_hash = {};
    my $record_id = '';
    my $parent_id = '';
    
    my @cols = split("\t");

    # convert to interbase coordinates:
    $cols[3]--;

    # change the + and - symbols in strand column to 0 and 1, respectively
    if ($cols[6] eq '+') {
        $cols[6] = 0;
    } elsif ($cols[6] eq '-') {
        $cols[6] = 1;
    } elsif ($cols[6] eq '.') {
        $cols[6] = 0;
        $logger->warn("converting feature with strand='.' to BSML feature with complement=false on line $linenum");
    } else {
        die("unknown value ($cols[6]) in strand column on line $linenum: expected '+', '-', or '.'.");
    }

    $record->{'_seqid'}     = uri_unescape($cols[0]);
    $record->{'_source'}    = uri_unescape($cols[1]);
    $record->{'_type'}      = uri_unescape($cols[2]);

    # URI-unescaping should not be necessary for these columns:
    $record->{'_start'}     = $cols[3];
    $record->{'_end'}       = $cols[4];
    $record->{'_score'}     = $cols[5];
    $record->{'_strand'}    = $cols[6];
    $record->{'_phase'}     = $cols[7];

    # swap start/end if necessary
    if ($record->{'_start'} > $record->{'_end'}) {
        $logger->warn("start > end on line $linenum, $cols[8]");
        $logger->logdie("start > end but strand is not '-' on line $linenum, $cols[8]") if ($record->{_strand} != 1);
        my $tmp = $record->{'_start'};
        $record->{'_start'} = $record->{'_end'};
        $record->{'_end'} = $tmp;
    }

    # apply $FEAT_TYPE_MAP to _type
    if (!defined($FEAT_TYPE_MAP->{$record->{'_type'}})) {
        $logger->error("unexpected feature type '$record->{_type}'");
    } else {
        $record->{'_type'} = $FEAT_TYPE_MAP->{$record->{'_type'}};
    }

    my @attribs = split(";", $cols[8]);
    foreach my $attrib(@attribs) {
        my ($type, $val) = split("=", $attrib);
        $type = uri_unescape($type);

        if (defined($unescaped_comma_atts->{$type})) {
            # treat $val as a single value, but other parts of the code assume that all attribute values are in lists
            $record->{$type}=[uri_unescape($val)];
        } else {
            my @vals = map { uri_unescape($_) } split(",", $val);
            $record->{$type}=\@vals;
        }
    }

    if (!defined($record->{'ID'})) {
        $logger->debug("no ID defined for record $line");
    } else {
        $record->{'ID'} = $record->{'ID'}->[0];
    }

    my $parents = [];

    if (!defined($record->{'Parent'})) {
        if (defined($record->{'Derives_from'})) {
            push(@$parents, @{$record->{'Derives_from'}});
        }
        $logger->debug("no Parent defined for record $line");
    } else {
        push(@$parents, @{$record->{'Parent'}});
    }

    foreach my $parent(@$parents) {
        # store record's hash reference as a child of Parent
        if (!defined($NODES->{$parent}->{'children'})) {
            $NODES->{$parent}->{'children'}->{$record->{'_type'}} = [];
        }
        push (@{$NODES->{$parent}->{'children'}->{$record->{'_type'}}},$record);
    }
    
    return $record;
}

# recursively traverse a node and all its children
# and process each record to create a feature hash
#
# $node - 
# $features - 
#
sub process_node {
    my ($node, $features) = @_;

    # hash of keys to ignore when extracting a record hash from a node hash
    my $ignore_keys = {
                        'children' => 1,
                      };
    
    # build a record hash from values stored in the node hash
    my $record;
    foreach my $key(keys %{$node}) {
        if (!$ignore_keys->{$key}) {
            $record->{$key} = $node->{$key};
        }
    }
    
    # process the record hash
    process_record($record, $features);
    
    # process the node's children
    foreach my $child_type(keys %{$node->{'children'}}) {
        foreach my $child_record(@{$node->{'children'}->{$child_type}}) {
                process_node($child_record, $features); 
        }
    }
    
    return; 
}

# process the records and store them in the features hash
sub process_record {
    my ($record, $features) = @_;
    my $cds_type = $FEAT_TYPE_MAP->{'CDS'};
    
    if ($record->{'_type'} eq $cds_type) {
        # CDS records can span lines and must be merged into one CDS feature
        if (!$record->{'ID'}) {
            if ($options->{'use_cds_parent_ids'}) {
                if ($record->{'Parent'}) {
                    my $trans_id = &get_parent_id_by_type({'parent' => $record->{'Parent'}}, $features, $MRNA_TYPES);
                    if (!defined($trans_id)) {
                        print Dumper $record;
                        die "no mRNA parent id found for CDS feature: try running with --insert_missing_mrnas";
                    }
                    $record->{'ID'} = $trans_id . "-CDS";
                } else {
                    die "CDS feature lacks ID and also has no Parent ID for --use_cds_parent_ids!";
                }
            } else {
                die "CDS feature lacks ID -> bad form!";
            }
        }
        if ($features->{$cds_type}->{$record->{'ID'}}) {
            if ($features->{$cds_type}->{$record->{'ID'}}->{'startpos'} > $record->{'_start'}) {
                $features->{$cds_type}->{$record->{'ID'}}->{'startpos'} = $record->{'_start'};
            }
            if ($features->{$cds_type}->{$record->{'ID'}}->{'endpos'} < $record->{'_end'}) {
                $features->{$cds_type}->{$record->{'ID'}}->{'endpos'} = $record->{'_end'};
            }
        } else {
            $features->{$cds_type}->{$record->{'ID'}} = {
                'complement'    => $record->{'_strand'},
                'startpos'      => $record->{'_start'},
                'endpos'        => $record->{'_end'},
                'parent'        => $record->{'Parent'},
                'type'          => $record->{'_type'},
                '_record'       => $record,
            };
        }
        
    } else {
        #handle all other feature types    

        my $id;
        if (!$record->{'ID'}) {
            $id = getTempId();
        } else {
            $id = $record->{'ID'};
        }
        
        my $title = undef;
        if ($record->{'Name'}) {
            my $name = shift(@{$record->{'Name'}});
            $title = $name if ($name !~ /^temp_id/);
        }

        $features->{$record->{'_type'}}->{$id} = {  
                    'complement'    => $record->{'_strand'},
                    'startpos'      => $record->{'_start'},
                    'endpos'        => $record->{'_end'},
                    'title'         => $title,
                    'parent'        => $record->{'Parent'} || $record->{'Derives_from'},
                    'type'          => $record->{'_type'},
                    '_record'       => $record,
                };
    }
    
}

# find the parent of a feature, checking only those parents in a specified
# set of parent types.  returns undef if a unique parent of one of the specified 
# types cannot be found
#
# $feat - feature whose parent to find
# $features - all features
# $types - hashref whose keys are the types of parents to check
#
sub get_parent_id_by_type {
    my($feat, $features, $types) = @_;
    my $parent_list = $feat->{'parent'};
    my $matches = [];

    foreach my $type (keys %$types) {
        my $type_list = $features->{$type};
        foreach my $parent_id (@$parent_list) {
            if (defined($type_list->{$parent_id})) {
                push(@$matches, $parent_id);
            }
        }
    }

    my $nm = scalar(@$matches);
    if ($nm != 1)  {
        $logger->warn("multiple parents found for feature $feat->{title}: " . join(',', @$parent_list)) if ($nm > 1);
        return undef;
    }

    return $matches->[0];
}

# convert a gene feature hash into BSML
sub gene_feature_hash_to_bsml {
    my ($nodes, $features, $gff_seqs, $peptides, $dna_seqs, $seq_id) = @_;
    
    my $id_hash = {};
    my %id_counts;

    $logger->debug("gene_feature_hash_to_bsml");
    foreach my $type(keys(%{$features})) {
        my $tc = scalar(keys(%{$features->{$type}}));
        $id_counts{$type} = $tc;
        $logger->debug("$type: $tc");
    }
    $idcreator->set_pool_size(%id_counts);
    
    my $seq;
    my $bsml_seq_id = $GFF2BSML_SEQID_MAP->{$seq_id};
    
    # create a sequence stub for the seq_id if it doesn't exist yet
    if (!defined($bsml_seq_id)) {

        # determine length of sequence from feature in GFF, if present
        my $gff_seq = $nodes->{$seq_id};
        my $gff_seqlen = undef;
        my $gff_seq_record = undef;

        if (defined($gff_seq)) {
            my $records = $gff_seq->{'records'};
            my $nr = scalar(@$records);
            if ($nr != 1) {
                $logger->warn("$nr GFF record(s) found for sequence with id = $seq_id");
            } else {
                # interbase coords in effect here
                $gff_seq_record = $records->[0];
                $gff_seqlen = $gff_seq_record->{'_end'} - $gff_seq_record->{'_start'};
            }
        }

        my $dseq = undef;
        my $dseq_file = undef;

        if (defined($dna_seqs) || defined($gff_seqs)) {
            $dseq = $dna_seqs->{$seq_id};
            $dseq_file = $options->{'dna_fasta'};
            my $dseq_i = defined($gff_seqs) ? $gff_seqs->{$seq_id} : undef;

            # look in both --dna_fasta and the file containing the sequences from the GFF itself
            if (defined($dseq) && defined($dseq_i)) {
                $logger->logdie("ambiguous sequence id: '$seq_id' was found in both the input GFF and --dna_fasta");
            }

            # use the sequence from the GFF if --peptide_fasta does not define one
            if (!defined($dseq)) {
                $dseq = $dseq_i;
                $dseq_file = $gff_seq_file;
            }

            die "couldn't find DNA sequence for $seq_id" if (!defined($dseq));
            $gff_seqlen = $dseq->{'seqlen'};
        }
        if (!defined($gff_seqlen)) {
            if (defined($dna_seqs)) {
                $logger->warn("unable to determine length of Sequence $seq_id: not found in GFF or --dna_fasta file");
            } else {
                $logger->warn("unable to determine length of Sequence $seq_id: feature not found in GFF and --dna_fasta not specified");
            }
        }

        my $bsml_id = $idcreator->next_id('project' => $options->{'project'}, 'type' => $options->{'default_seq_class'});
        $GFF2BSML_SEQID_MAP->{$seq_id} = $bsml_id;
        $seq = $doc->createAndAddSequence( $bsml_id, $seq_id, $gff_seqlen, 'dna', $options->{'default_seq_class'} );
        my $genome_link = $doc->createAndAddLink($seq, 'genome', '#' . $genome_id);

        my $default_atts = $options->{'default_seq_attributes'};
        if (defined($default_atts) && ($default_atts =~ /:/)) {
            foreach my $def_att (split(/,/, $default_atts)) {
                my($key, $val) = split(/:/, $def_att);
                if (defined($key) && defined($val)) {
                    $doc->createAndAddBsmlAttribute($seq, $key, $val);
                }                    
            }
        }

        if (defined($dseq)) {
            if (!$options->{'dna_no_seqdata'}) {
                my $path = File::Spec->rel2abs($dseq_file);
                my $seq_data_import_elem = $doc->createAndAddSeqDataImport($seq, 'fasta', $path, undef, $dseq->{'bsml_identifier'});
            }
            $doc->createAndAddBsmlAttribute($seq, 'defline', $dseq->{'defline'});
        }

        # handle cross-references
        foreach my $att (keys %$seq_xref_mappings) {
            my $mappings = $seq_xref_mappings->{$att};
            my $att_val;
            if (defined($gff_seq_record)) {
                $att_val = $gff_seq_record->{$att};
                unshift(@$att_val, $gff_seq_record->{'title'}) if (($att eq 'Name') && (defined($gff_seq_record->{'title'})) && ($gff_seq_record->{'title'} =~ /\S+/));
            } elsif ($att eq 'ID') {
                $att_val = $seq_id;
            }
            next unless (defined($att_val));
            my $att_vals = ((ref $att_val) eq 'ARRAY') ? $att_val : [$att_val];

            foreach my $mapping (@$mappings) {
                foreach my $av (@$att_vals) {
                    my $xref_elem = $doc->createAndAddCrossReference('parent'          => $seq,
                                                                     'id'              => $doc->{'xrefctr'}++,
                                                                     'database'        => $mapping->{'db_name'},
                                                                     'identifier'      => $av,
                                                                     'identifier-type' => $mapping->{'id_type'},
                                                                     );
                }
            }
        }

        # TODO - need a way to specify which features get linked to an analysis
#        $seq->addBsmlLink('analysis', '#EVM', 'input_of');
    } else {
        $seq = $doc->returnBsmlSequenceByIDR($bsml_seq_id);
    }

    # each of the feature types listed in $MRNA_TYPES can be used as a transcript-level feature
    my @transcript_ids = ();
    foreach my $ttype (keys %$MRNA_TYPES) {
        my @mrna_ids = map {{'id' => $_, 'type' => $ttype}} keys(%{$features->{$ttype}});
        push(@transcript_ids, @mrna_ids);
    }

    # create at most one feature table per sequence
    my $feat_table = $seq2feat_table->{$seq};
    if (!defined($feat_table)) {
        $feat_table = $seq2feat_table->{$seq} = $doc->createAndAddFeatureTable($seq);
    }

    # $feat_groups indexed by transcript id
    my $feat_groups = {};

    foreach my $trans_id(@transcript_ids) {
        my $t_id = $idcreator->next_id( 
                                        'project' => $options->{'project'}, 
                                        'type' => $trans_id->{'type'},
                                        );

        if (defined($feat_groups->{$trans_id->{'id'}})) {
            die "duplicate transcript id ($trans_id->{id}) encountered";
        } else {
            $id_hash->{$trans_id->{'id'}} = $t_id;
            $feat_groups->{$trans_id->{'id'}} = $doc->createAndAddFeatureGroup($seq, '', $t_id);
            
            # create reference for parent of transcript also
            my $trans = $features->{$trans_id->{'type'}}->{$trans_id->{'id'}};
            foreach my $parent_id (@{$trans->{'parent'}}) {
                $feat_groups->{$parent_id} = $feat_groups->{$trans_id->{'id'}};
            }
        }
    }

    my $get_singleton_feat = sub {
        my($type, $return_arrays) = @_;
        my @feats = values %{$features->{$type}};
        my @ids = keys %{$features->{$type}};
        my $nf = scalar(@feats);
        # not using wantarray because the return value is always an array/pair
        if ($return_arrays) {
            return (\@ids, \@feats);
        }
        if ($nf != 1) {
            print Dumper $features;
            die "unable to find unique $type feature associated with $seq_id ($nf $type feature(s) found)";
        }
        return ($ids[0], $feats[0]);
    };

    my $cds_type = $FEAT_TYPE_MAP->{'CDS'};
    my $pep_type = $FEAT_TYPE_MAP->{'polypeptide'};
    my $added_pep = 0;
    if ($options->{'insert_polypeptides'} && (!defined($id_counts{$pep_type}))) {
        my($gene_id, $gene) = &$get_singleton_feat('gene');
        my($trans_ids, $transs) = &$get_singleton_feat('transcript', 1);

        # check whether this is a protein-coding gene
        my $is_protein_coding = 1;
        my $nc_mrna_type = undef;
        my $nc_mrna_count = 0;

        if (scalar(@$trans_ids) == 0) { 
            foreach my $mrna_type (keys %$NONCODING_RNA_TYPES) {
                my($ids, $ts) = &$get_singleton_feat($mrna_type, 1);
                if (scalar(@$ids) == 1) {
                    $nc_mrna_type = $mrna_type;
                    ++$nc_mrna_count;
                }
            }
            if ($nc_mrna_count == 1) {
                $is_protein_coding = 0;
            } else {
                die "gene $gene_id may be noncoding but has $nc_mrna_count noncoding mRNA feature(s)";
            }
        } elsif (scalar(@$trans_ids) > 1) {
                die "gene $gene_id has multiple mRNA feature(s)";
        }

        my($trans_id, $trans) = ($trans_ids->[0], $transs->[0]);

        if ($is_protein_coding) {
            my($cds_ids, $cdss) = &$get_singleton_feat($cds_type, 1);
            if (scalar(@$cds_ids) > 1) {
                $logger->logdie("gene $gene_id has multiple CDS feature(s) after merging");
            } 
            # not protein coding after all (or an error)
            elsif (scalar(@$cds_ids) == 0) {
                my $msg = "gene $gene_id has no CDS feature and it does not appear to be an RNA gene";
                if ($options->{'allow_genes_with_no_cds'}) { 
                    $logger->warn($msg); 
                    $is_protein_coding = 0;
                }
                else { 
                    $logger->logdie($msg); 
                }
            }

            if ($is_protein_coding) {
                my($cds_id, $cds) = ($cds_ids->[0], $cdss->[0]);
                my $pep_id = undef;
                my $id_type = $options->{'inserted_polypeptide_id_type'};
                if ($id_type eq 'gene') {
                    $pep_id = $gene_id . '-Protein';
                } elsif ($id_type eq $cds_type) {
                    $pep_id = $cds_id . '-Protein';
                } else { # mrna/transcript
                    $pep_id = $trans_id . '-Protein';
                }
                
                $features->{$pep_type}->{$pep_id} = {
                    'parent' => [ $cds_id ],
                    'type' => $pep_type,
                    'startpos' => $cds->{'startpos'},
                    'endpos' => $cds->{'endpos'},
                    'complement' => $cds->{'complement'},
                    # TODO - propagation of cross-references and attributes to automatically-inserted features should be optional/explicit:
                    '_record' => $cds->{'_record'},
                };
                $added_pep = 1;
            }
        }
    }

    # track BSML features for attribute/xref cloning
    my $id2bsmlfeat = {};

    foreach my $type(keys(%{$features})) {
        foreach my $feat_id(keys(%{$features->{$type}})) {
                my $id;
                if (! defined($id_hash->{$feat_id})) {
                    $id = $idcreator->next_id( 
                                               'project' => $options->{'project'}, 
                                               'type' => $type 
                                               );
                } else {
                    $id = $id_hash->{$feat_id};
                }

                my $feat_title = ($feat_id =~ /temp_id/) ? undef : $feat_id;
                my $feat = $doc->createAndAddFeature( $feat_table, $id, $feat_title, $type);
                $id2bsmlfeat->{$feat_id} = $feat;

                # link polypeptides to Sequences in --peptide_fasta or the input GFF file
                if (($type eq 'polypeptide') && (defined($peptides) || defined($gff_seqs))) {
                    my $fseq = $peptides->{$feat_id};
                    my $fseq_file = $options->{'peptide_fasta'};
                    my $fseq_i = defined($gff_seqs) ? $gff_seqs->{$feat_id} : undef;

                    # look in both --peptide_fasta and the file containing the sequences from the GFF itself
                    if (defined($fseq) && defined($fseq_i)) {
                        $logger->logdie("ambiguous sequence id: '$feat_id' was found in both the input GFF and --peptide_fasta");
                    }

                    # use the sequence from the GFF if --peptide_fasta does not define one
                    if (!defined($fseq)) {
                        $fseq = $fseq_i;
                        $fseq_file = $gff_seq_file;
                    }

                    if (!defined($fseq)) {
                        if ($added_pep) {
                            $logger->warn("couldn't find peptide sequence for '$feat_id'");
                        } else {
                            die "couldn't find peptide sequence for '$feat_id'";
                        }
                    } else {
                        my $pep_seq = $doc->createAndAddSequence($id . "_seq", undef, $fseq->{'seqlen'}, 'aa', 'polypeptide');
                        my $seq_genome_link = $doc->createAndAddLink($pep_seq, 'genome', '#' . $genome_id);
                        my $feat_seq_link = $doc->createAndAddLink($feat, 'sequence', '#'.$id.'_seq');
                        if (!$options->{'peptide_no_seqdata'}) {
                            my $path = File::Spec->rel2abs($fseq_file);
                            my $seq_data_import_elem = $doc->createAndAddSeqDataImport($pep_seq, 'fasta', $path, undef, $fseq->{'bsml_identifier'});
                        }
                        $doc->createAndAddBsmlAttribute($pep_seq, 'defline', $fseq->{'defline'});
                    }
                }

                # TODO - need a way to specify which features get linked to an analysis
#                $feat->addBsmlLink('analysis', '#EVM', 'computed_by');

                $feat->addBsmlIntervalLoc(
                                          $features->{$type}->{$feat_id}->{'startpos'},
                                          $features->{$type}->{$feat_id}->{'endpos'},
                                          $features->{$type}->{$feat_id}->{'complement'},
                                          );

                # add attributes
                foreach my $att (keys %$feat_att_mappings) {
                    my $feature = $features->{$type}->{$feat_id};
                    my $mappings = $feat_att_mappings->{$att};
                    my $att_val = $feature->{'_record'}->{$att};
                    unshift(@$att_val, $feature->{'title'}) if (($att eq 'Name') && (defined($feature->{'title'})) && ($feature->{'title'} =~ /\S+/));
                    
                    if (defined($att_val)) {
                        my $att_vals = ((ref $att_val) eq 'ARRAY') ? $att_val : [$att_val];

                        foreach my $mapping (@$mappings) {
                            foreach my $av (@$att_vals) {
                                $doc->createAndAddBsmlAttribute($feat, $mapping->{'bsml_att'}, $av);
                            }
                        }
                    }
                }

                # add cross-references
                foreach my $att (keys %$feat_xref_mappings) {
                    my $feature = $features->{$type}->{$feat_id};
                    my $mappings = $feat_xref_mappings->{$att};
                    my $att_val = $feature->{'_record'}->{$att};
                    unshift(@$att_val, $feature->{'title'}) if (($att eq 'Name') && (defined($feature->{'title'})) && ($feature->{'title'} =~ /\S+/));
                    
                    if (defined($att_val)) {
                        my $att_vals = ((ref $att_val) eq 'ARRAY') ? $att_val : [$att_val];

                        foreach my $mapping (@$mappings) {
                            foreach my $av (@$att_vals) {
                                my $xref_elem = $doc->createAndAddCrossReference('parent'          => $feat,
                                                                                 'id'              => $doc->{'xrefctr'}++,
                                                                                 'database'        => $mapping->{'db_name'},
                                                                                 'identifier'      => $av,
                                                                                 'identifier-type' => $mapping->{'id_type'},
                                                                                 );
                            }
                        }
                    }
                }

                # determine which feat_group this feature belongs with
                my $feat_group = undef;
                if (defined($feat_groups->{$feat_id})) { 
                    $feat_group = $feat_groups->{$feat_id};
                } else {
                    # polypeptides _should_ be linked through CDS
                    my $cds_id = undef;
                    my $cds_type = undef;
                    if ($type eq $pep_type) {
                        $cds_type = $FEAT_TYPE_MAP->{'CDS'};
                        $cds_id = &get_parent_id_by_type($features->{$type}->{$feat_id}, $features, { $cds_type => 1 });
                    }

                    my $trans_id = undef;
                    if (defined($cds_id)) {
                        $trans_id = &get_parent_id_by_type($features->{$cds_type}->{$cds_id}, $features, $MRNA_TYPES);
                    } else {
                        $trans_id = &get_parent_id_by_type($features->{$type}->{$feat_id}, $features, $MRNA_TYPES);
                    }

                    if (!defined($trans_id)) {
                        print Dumper $features;
                        die "no mRNA parent id found for $type feature $feat_id";
                    }
                    $feat_group = $feat_groups->{$trans_id};
                    if (!defined($feat_group)) {
                        print Dumper $features;
                        die "no feat_group found for transcript=$trans_id";
                    }
                }
                
                $feat_group->addBsmlFeatureGroupMember( $id, $type );
            }
    }

    # clone BSML attributes specified by --clone_feat_attribute_mappings, --clone_feat_xref_mappings
    my $all_clone_mappings = {};
    foreach my $key (keys %$clone_feat_att_mappings, keys %$clone_feat_xref_mappings) {
        next if defined($all_clone_mappings->{$key});
        my $fam = $clone_feat_att_mappings->{$key};
        my $fxm = $clone_feat_xref_mappings->{$key};
        map { $_->{'mapping_type'} = 'attribute'; } @$fam;
        map { $_->{'mapping_type'} = 'xref'; } @$fxm;
        $all_clone_mappings->{$key} = [@$fam, @$fxm];
    }

    foreach my $st (keys %$all_clone_mappings) {
        my $mappings = $all_clone_mappings->{$st};
        # get features of type $st
        my @src_feat_ids = keys %{$features->{$st}};
        $logger->debug("got src feats of type $st: @src_feat_ids");

        foreach my $cfm (@$mappings) {
            my($att, $xref_db, $tt, $mt) = map {$cfm->{$_}} ('bsml_att', 'db_name', 'target_type', 'mapping_type');
            # get features of type $tt
            my @tgt_feat_ids = keys %{$features->{$tt}};
            $logger->debug("got tgt feats of type $tt: @tgt_feat_ids");
            # clone
            foreach my $sfid (@src_feat_ids) {
                foreach my $tfid (@tgt_feat_ids) {
                    my $sf = $id2bsmlfeat->{$sfid} || $logger->logdie("no feature found for id=$sfid");
                    my $tf = $id2bsmlfeat->{$tfid} || $logger->logdie("no feature found for id=$tfid");
                    
                    # clone attribute
                    if ($mt eq 'attribute') {
                        $logger->debug("cloning $att from $st to $tt");
                        my $tvh = {}; # track target atts to avoid duplication
                        my ($sval, $tval) = map { $_->returnBsmlAttr($att) || []; } ($sf, $tf);
                        map {$tvh->{$_} = 1;} @$tval;
                        foreach my $sv (@$sval) {
                            if (!defined($tvh->{$sv})) {
                                $doc->createAndAddBsmlAttribute($tf, $att, $sv);
                                $logger->debug("cloning $att value of $sv from $sfid ($st) to $tfid ($tt)");
                            }
                        }
                    } 
                    # clone cross-reference
                    elsif ($mt eq 'xref') {
                        $logger->debug("cloning $xref_db xref from $st to $tt");
                        my $trefh = {}; # track target xrefs to avoid duplication
                        my $trefs = $tf->returnBsmlCrossReferenceListR();
                        foreach my $tr (@$trefs) {
                            my($tdb,$tid_type,$tid) = map { $tr->returnattr($_) || ""; } ('database', 'identifier-type', 'identifier');
                            $trefh->{join("\0",$tdb,$tid_type,$tid)} = 1;
                        }

                        # TODO - make this more efficient
                        my $srefs = $sf->returnBsmlCrossReferenceListR();

                        foreach my $sr (@$srefs) {
                            my($sdb,$sid_type,$sid) = map { $sr->returnattr($_); } ('database', 'identifier-type', 'identifier');
                            next unless ($sdb eq $xref_db);

                            # add it if the target doesn't have it already
                            if (!defined($trefh->{join("\0", $sdb,$sid_type,$sid)})) {
                                my $cind = $tf->addBsmlCrossReference();
                                my $cr = $tf->returnBsmlCrossReferenceR($cind);
                                $cr->addattr('database', $sdb);
                                $cr->addattr('identifier-type', $sid_type);
                                $cr->addattr('identifier', $sid);
                                $logger->debug("adding cross ref from target feature $tfid to db=$sdb idtype=$sid_type, id=$sid");
                            } else {
                                $logger->debug("target feature $tfid already has cross ref with db=$sdb idtype=$sid_type, id=$sid from source feature $sfid");
                            }
                        }

                    } else {
                        $logger->logdie("attempt to clone unknown BSML component '$mt'");
                    }
                }
            }
        }
    }

    # TODO - clone BSML cross-references specified by --clone_feat_xref_mappings
    
}

# traverses an array of node references 
# and returns all nodes of the specified type
sub fetch_node_type {
    my ($type, $nodes_ref, $found_nodes) = @_;
    
    foreach my $node (@{$nodes_ref}) {
        if ($node->{'_type'} eq $type) {
            push (@{$found_nodes}, $node);
        }
        foreach my $key (keys %{$node->{'children'}}) {
            fetch_node_type($type, $node->{'children'}->{$key}, $found_nodes);
        }
    }
    
    return;
}

# read deflines and sequence lengths from a multi-FASTA file
#
# $file - a multi-FASTA sequence file
# $id_regex - regular expression to extract a sequence id from each defline of $file
# $id_fn - function to apply to the id parsed by $id_regex before using it as a key in the result hash
#
# returns a hashref mapping sequence id to { 'seqlen' => sequence length, 'defline' => sequence defline }
#
sub read_sequence_lengths {
    my($file, $id_regex, $id_fn) = @_;
    my $result = {};
    my $fh = FileHandle->new();
    my $bsml_idents = {};

    my $process_seq = sub {
        my($defline, $seq, $lnum) = @_;
        my($id) = ($defline =~ /$id_regex/);
        die "unable to parse id from line $lnum of $file using regex '$id_regex': $defline" if (!defined($id));
        $seq =~ s/\s+//g;
        $id = &$id_fn($id) if (defined($id_fn));
        die "sequence id '$id' is not unique in $file" if (defined($result->{$id}));
        # this must mirror the hard-coded regex in BSML::BsmlReader::parse_multi_fasta:
        my($bsml_ident) = ($defline =~ /^>([^\s]+)/);
        die "unable to parse BSML identifier from the following defline in $file: $defline" if (!defined($bsml_ident));
        die "BSML identifier '$bsml_ident' is not unique in $file" if (defined($bsml_idents->{$bsml_ident}));
        $bsml_idents->{$bsml_ident} = 1;
        $defline =~ s/^>//;
        $result->{$id} = { 'id' => $id, 'seqlen' => length($seq), 'defline' => $defline, 'bsml_identifier' => $bsml_ident };
    };

    my $defline = undef;
    my $defline_lnum = undef;
    my $seq = undef;
    my $lnum = 0;
    my $openCmd = $file;
    if ($file =~ /\.gz$/) {
        $openCmd = "zcat $file |";
        $logger->warn("FASTA file $file is compressed, which may cause bsml2chado to fail to load sequences from it.");
    }

    $fh->open($openCmd) || die "unable to open $openCmd for file $file";
    while (my $line = <$fh>) {
        chomp($line);
        ++$lnum;
        if ($line =~ /^>/) {
            &$process_seq($defline, $seq, $defline_lnum) if (defined($defline));
            $defline = $line;
            $defline_lnum = $lnum;
            $seq = "";
        } else {
            $seq .= $line;
        }
    }
    &$process_seq($defline, $seq, $defline_lnum) if (defined($defline));
    $fh->close();

    return $result;
}

sub getTempId {
    return "temp_id_".$GLOBAL_ID_COUNTER++; 
}
