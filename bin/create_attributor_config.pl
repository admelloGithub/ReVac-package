#!/usr/local/bin/perl -w
use strict;
use Getopt::Std;

getopts("h:H:R:T:L:F:");
our ($opt_h,$opt_H,$opt_R,$opt_T,$opt_L,$opt_F);

if ($opt_h) {
	print "This script generates a config file for the attributor component.\n
			Options are:\n
			-H for Hmmpfam3 htab list\n
			-R for RAPSearch2 m8 formatted list\n
			-T for TMHMM output list\n
			-L for LipoP bsml list\n
			-F for fasta input used in generating above outputs\n";			
}

if (defined $opt_H && defined $opt_R && defined $opt_T && defined $opt_L && defined $opt_F) {

print "
---
#Config file for IGS Prokaryotic Functional Annotation 
#According to guidlines here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3111993/pdf/sigs.1223234.pdf
 
general:
   # All proteins will start out with this product name, before evidence
   #  is evaluated
   default_product_name: hypothetical protein
   
   # If Yes, each annotation attribute is allowed from a different source.  For
   #  example, if there is a high-quality HMM hit that provides a gene product
   #  name and GO terms but lacks an EC number, another source (such as a BLAST
   #  hit) may provide the EC number.
   allow_attributes_from_multiple_sources: No

   # This is purely for development purposes and should usually best set to 0.  If
   #  you set this to any other integer, it will limit the number of polypeptides
   #  for which results are considered.  This allows for more rapid testing of
   #  the parser on larger datasets.
   debugging_polypeptide_limit: 0

indexes:
   coding_hmm_lib: /usr/local/projects/PNTHI/databases/attributor/dbs/coding_hmm_lib.sqlite3
   uniref100: /local/projects/uniref/uniprot_trembl.sqlite3
   uniprot_sprot: /usr/local/projects/PNTHI/databases/attributor/dbs/uniprot_sprot.20171025.sqlite3

input:
   # These are the files on which annotation will be applied.  At a minimum, the 'fasta'
   #  option must be defined.  GFF3 output cannot be specified unless 'gff3' input is
   #  also provided.
   # These test data are for SRS15430
   polypeptide_fasta: $opt_F 
   #gff3: /local/projects/PFDA1/snadendla/batch4/AMERTCC_69/genemark_es/AMERTCC69_genemark_es.gff3/genemark.gff3

order:
   - coding_hmm_lib__equivalog
   - rapsearch2__uniref100__trusted_full_full
   - coding_hmm_lib__equivalog_domain
   - rapsearch2__uniref100__trusted_partial_full
   - coding_hmm_lib__subfamily
   - coding_hmm_lib__superfamily
   - coding_hmm_lib__subfamily_domain
   - coding_hmm_lib__domain
   - coding_hmm_lib__pfam
   - rapsearch2__uniref100__trusted_full_partial
   - tmhmm
   - coding_hmm_lib__hypothetical_equivalog
#   - lipoprotein_motif

evidence:
   - label: coding_hmm_lib__equivalog
     type: HMMer3_htab
     path: $opt_H
     class: equivalog
     index: coding_hmm_lib

   - label: coding_hmm_lib__equivalog_domain
     type: HMMer3_htab
     path: $opt_H
     class: equivalog_domain
     index: coding_hmm_lib

   - label: coding_hmm_lib__subfamily
     type: HMMer3_htab
     path: $opt_H
     class: subfamily
     index: coding_hmm_lib
     append_text: family protein

   - label: coding_hmm_lib__superfamily
     type: HMMer3_htab
     path: $opt_H
     class: superfamily
     index: coding_hmm_lib
     append_text: family protein

   - label: coding_hmm_lib__subfamily_domain
     type: HMMer3_htab
     path: $opt_H
     class: subfamily_domain
     index: coding_hmm_lib
     append_text: domain protein

   - label: coding_hmm_lib__domain
     type: HMMer3_htab
     path: $opt_H
     class: domain
     index: coding_hmm_lib
     append_text: domain protein

   - label: coding_hmm_lib__pfam
     type: HMMer3_htab
     path: $opt_H
     class: pfam
     index: coding_hmm_lib
     append_text: family protein

   - label: coding_hmm_lib__hypothetical_equivalog
     type: HMMer3_htab
     path: $opt_H
     class: hypoth_equivalog
     index: coding_hmm_lib

   - label: rapsearch2__uniref100__trusted_full_full
     type: RAPSearch2_m8
     path: $opt_R
     class: trusted
     index: uniref100
     query_cov: 80%
     match_cov: 80%
     percent_identity_cutoff: 50%

   - label: rapsearch2__uniref100__trusted_partial_full
     type: RAPSearch2_m8
     path: $opt_R
     class: trusted
     index: uniref100
     match_cov: 80%
     percent_identity_cutoff: 50%
     append_text: domain protein

   - label: rapsearch2__uniref100__trusted_full_partial
     type: RAPSearch2_m8
     path: $opt_R
     class: trusted
     index: uniref100
     query_cov: 80%
     percent_identity_cutoff: 50%
     append_text: domain protein

   - label: rapsearch2__uniref100__all_full_full
     type: RAPSearch2_m8
     path: $opt_R
     index: uniref100
     query_cov: 80%
     match_cov: 80%
     percent_identity_cutoff: 50%
     prepend_text: putative

   - label: rapsearch2__uniref100__all_partial_full
     type: RAPSearch2_m8
     path: $opt_R
     index: uniref100
     match_cov: 80%
     percent_identity_cutoff: 50%
     prepend_text: putative
     append_text: domain protein

   - label: rapsearch2__uniref100__all_full_partial
     type: RAPSearch2_m8
     path: $opt_R
     index: uniref100
     query_cov: 80%
     percent_identity_cutoff: 50%
     prepend_text: putative
     append_text: domain protein

   - label: tmhmm
     type: TMHMM
     # this is the product name that will be assigned given a positive TMHMM match
     product_name: putative integral membrane protein
     # minimum required predicted helical spans across the membrane required to apply evidence
     min_helical_spans: 5
     #path: /usr/local/projects/dacc/output_repository/tmhmm/14365_default/tmhmm.raw.list
     #path: /usr/local/projects/aengine//output_repository/tmhmm/12808815406_default/tmhmm.raw.list
     path: $opt_T

#   - label: lipoprotein_motif
#     type: lipoprotein_motif_bsml
#     path: $opt_L
     # This is the product name that will be assigned given a positive match
#     product_name: Putative lipoprotein
";

} else {
	
	die "\nNeed all arguments -HRTLF\nUse -h for details.\n";

}
#EOF
