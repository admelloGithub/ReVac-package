#!/usr/bin/env python3

"""
This is to assign functional annotation to gene models using a wide variety of evidence.  Meant to
support unlimited hierarchy rules, this utility relies on a user-created configuration file.
Current limitations:
- Rules are evaluated on a per-gene basis, so no rules are currently possible that
  would need to consider annotations of other genes in the set.
- IGS 'PFunc' hierarchy
http://imgur.com/1odYcT5
"""

import argparse
import cProfile
import math
import os
import re

import biocode.annotation
import biocode.gff
import biocode.utils
import biocode.things

import sqlite3
import sys
import xml.etree.ElementTree as etree
import yaml

outfile = sys.argv[-1]
outdir = outfile.split('/')
d = '/'
out_dir = d.join(outdir[:-1])

def main():
    parser = argparse.ArgumentParser( description='Assigns functional annotation based on user-configurable evidence tiers')

    ## output file to be written
    parser.add_argument('-c', '--config_file', type=str, required=True, help='Configuration file for annotation' )
    parser.add_argument('-o', '--output_base', type=str, required=True, help='Base name/path of output files to be created' )
    parser.add_argument('-f', '--output_format', type=str, required=False, default='gff3', help='Desired output format' )
    args = parser.parse_args()

    sources_log_fh = open("{0}.sources.log".format(args.output_base), 'wt')

    configuration = yaml.load(open(args.config_file).read())
    check_configuration(configuration, args)
    evidence = parse_evidence_config(configuration)
    default_product_name = configuration['general']['default_product_name']

    # stores any active SQLite3 db connections
    db_conn = dict()

    # this is a dict of biothings.Polypeptide objects
    polypeptides = initialize_polypeptides(sources_log_fh, configuration['input']['polypeptide_fasta'], default_product_name)

    for label in configuration['order']:
        if label not in evidence:
            raise Exception("ERROR: There is a label '{0}' in the 'order' section of the conf file that isn't present in the 'evidence' section".format(label))

        if evidence[label]['type'] == 'HMMer3_htab':
            index_conn, ev_db_conn = get_or_create_db_connections(type_ev='hmm_ev', configuration=configuration,
                                         evidence=evidence, label=label, db_conn=db_conn, output_base=args.output_base)
            index_conn.isolation_level = None
            apply_hmm_evidence(polypeptides=polypeptides, ev_conn=ev_db_conn, config=configuration,
                               ev_config=evidence[label], label=label, index_conn=index_conn, log_fh=sources_log_fh)
                
        elif evidence[label]['type'] == 'RAPSearch2_m8':
            index_conn, ev_db_conn = get_or_create_db_connections(type_ev='blast_ev', configuration=configuration,
                                         evidence=evidence, label=label, db_conn=db_conn, output_base=args.output_base)
            index_conn.isolation_level = None
            apply_blast_evidence(polypeptides=polypeptides, ev_conn=ev_db_conn, config=configuration,
                                 ev_config=evidence[label], label=label, index_conn=index_conn, log_fh=sources_log_fh)

        elif evidence[label]['type'] == 'TMHMM':
            index_conn, ev_db_conn = get_or_create_db_connections(type_ev='tmhmm_ev', configuration=configuration,
                                         evidence=evidence, label=label, db_conn=db_conn, output_base=args.output_base)
            apply_tmhmm_evidence(polypeptides=polypeptides, ev_conn=ev_db_conn, config=configuration,
                                 ev_config=evidence[label], label=label, log_fh=sources_log_fh)
            
        elif evidence[label]['type'] == 'lipoprotein_motif_bsml':
            index_conn, ev_db_conn = get_or_create_db_connections(type_ev='lipoprotein_motif_ev', configuration=configuration,
                                         evidence=evidence, label=label, db_conn=db_conn, output_base=args.output_base)
            apply_lipoprotein_motif_evidence(polypeptides=polypeptides, ev_conn=ev_db_conn, config=configuration,
                                             ev_config=evidence[label], label=label, log_fh=sources_log_fh)

        else:
            raise Exception("ERROR: Unsupported evidence type '{0}' with label '{1}' in configuration file".format(evidence[label]['type'], label))

    # close all db connections
    for label in db_conn:
        db_conn[label].close()

    perform_final_checks(polypeptides=polypeptides, config=configuration, log_fh=sources_log_fh)

    # Write the output
    polyset = biocode.things.PolypeptideSet()
    polyset.load_from_dict(polypeptides)
    
    if args.output_format == 'fasta':
        polyset.write_fasta(path="{0}.faa".format(args.output_base))
    elif args.output_format == 'gff3':
        ## parse input GFF
        (assemblies, ref_features) = biocode.gff.get_gff3_features( configuration['input']['gff3'] )

        ## merge annotation with polypeptide collection
        biocode.gff.add_annotation(features=ref_features, polypeptide_set=polyset)

        ## print the new GFF
        biocode.gff.print_gff3_from_assemblies(assemblies=assemblies, ofh=open("{0}.gff3".format(args.output_base), 'wt'))

    
def already_indexed(path=None, index=None):
    curs = index.cursor()
    curs.execute("SELECT id FROM data_sources WHERE source_path = ?", (path, ))
    found = False

    for row in curs:
        found = True
        break
    
    curs.close()
    return found

def apply_blast_evidence(polypeptides=None, ev_conn=None, config=None, ev_config=None, label=None, index_conn=None, log_fh=None):
    """
    Uses BLAST (or similar) evidence to assign functional evidence to polypeptides.  Description of arguments:
    polypeptides: a dict of biothings.Polypeptides objects, keyed on ID
    ev_conn: SQLite3 connection to the parsed BLAST(ish) search evidence db for this set of searches
    config: The yaml object for the parsed annotation config file
    ev_config: The parsed evidence section for this label within the annotation config file
    label: Label for the evidence track entry within the annotation config file
    index_conn:  SQLite3 connection to the reference index for the database searched
    """
    print("DEBUG: Processing BLAST evidence tier label: {0}".format(label))
    default_product = config['general']['default_product_name']
    ev_curs = ev_conn.cursor()

    # Doing all the cursors here prevents doing it repeatedly within calling each accession.  Huge performance improvement.
    index_acc_curs = index_conn.cursor()
    index_acc_curs.execute("begin")
    index_ec_curs = index_conn.cursor()
    index_go_curs = index_conn.cursor()

    ev_qry = "SELECT sbj_id, align_len, perc_identity, eval, bit_score FROM blast_hit WHERE qry_id = ? ORDER BY eval ASC"

    print("DEBUG: Applying {1} results to {0} polypeptides".format(len(polypeptides), label))

    blast_class_limit = None
    
    if 'class' in ev_config:
        blast_class_limit = ev_config['class']

    prepend_text = None
    if 'prepend_text' in ev_config:
        prepend_text = ev_config['prepend_text']

    append_text = None
    if 'append_text' in ev_config:
        append_text = ev_config['append_text']
        
    if 'debugging_polypeptide_limit' in config['general']:
        DEBUG_LIMIT = config['general']['debugging_polypeptide_limit']

    # Are coverage cutoffs defined?
    query_cov_cutoff = None
    match_cov_cutoff = None
    if 'query_cov' in ev_config:
        query_cov_cutoff = int(ev_config['query_cov'].rstrip('%'))

    if 'match_cov' in ev_config:
        match_cov_cutoff = int(ev_config['match_cov'].rstrip('%'))

    percent_identity_cutoff = None
    if 'percent_identity_cutoff' in ev_config:
        percent_identity_cutoff = int(ev_config['percent_identity_cutoff'].rstrip('%'))

    for id in polypeptides:
        polypeptide = polypeptides[id]
        print("DEBUG: Parsing {0} evidence for polypeptide ID {1}, length: {2}".format(label, id, polypeptide.length))
        annot = polypeptide.annotation

        DEBUG_LIMIT = DEBUG_LIMIT - 1
        if DEBUG_LIMIT == 0:
            break

        if config['general']['allow_attributes_from_multiple_sources'] == 'Yes':
            raise Exception("ERROR: Support for the general:allow_attributes_from_multiple_sources=Yes setting not yet implemented")
        else:
            if annot.product_name != default_product: continue

            for ev_row in ev_curs.execute(ev_qry, (polypeptide.id,)):
                if query_cov_cutoff is not None:
                    perc_coverage = (ev_row[1] / polypeptide.length) * 100
                    if perc_coverage < query_cov_cutoff:
                        print("\tSkipping accession {0} because coverage {1}% doesn't meet cutoff {2}% requirement".format(
                            ev_row[0], perc_coverage, query_cov_cutoff))
                        continue

                if percent_identity_cutoff is not None:
                    if ev_row[2] < percent_identity_cutoff:
                        print("\tSkipping accession {0} because percent identity of {1}% doesn't meet cutoff {2}% requirement".format(
                            ev_row[0], ev_row[2], percent_identity_cutoff))
                        continue

                ## This is completely crap that we have to do this, and is a byproduct of the fact that the
                #   database metadata indexing tools create different table names 
                m = re.search('uniref', label, re.I)
                if m:
                    blast_annot = get_uniref_accession_info(conn=index_conn, accession=ev_row[0], config=config, acc_curs=index_acc_curs, ec_curs=index_ec_curs, go_curs=index_go_curs)
                else:
                    m = re.search('sprot', label, re.I)
                    if m:
                        blast_annot = get_sprot_accession_info(conn=index_conn, accession=ev_row[0], config=config)
                    else:
                        raise Exception("ERROR: Expected search index label '{0}' to contain either 'uniref' or 'sprot' in order to know which metadatadb to use.".format(label))

                # is there a class limitation?
                if blast_class_limit is not None:
                    if blast_class_limit == 'trusted':
                        if 'is_characterized' in blast_annot.other_attributes:
                            if blast_annot.other_attributes['is_characterized'] != 1:
                                print("\tSkipping accession {0} because it is not characterized".format(ev_row[0]))
                                continue
                            else:
                                print("\tAccepting accession {0} because it is characterized".format(ev_row[0]))
                                pass
                    else:
                        raise Exception("ERROR: Unrecognized value ('{0}') for class in config file".format(blast_class_limit))
                    
                if match_cov_cutoff is not None:
                    if 'ref_len' not in blast_annot.other_attributes:
                        print("\tSkipping accession {0} because length wasn't found".format(ev_row[0]))
                        continue
                    
                    match_coverage = (ev_row[1] / blast_annot.other_attributes['ref_len'])*100
                    if match_coverage < match_cov_cutoff:
                        print("\tSkipping accession {0} because match coverage {1}% doesn't meet cutoff {2}% requirement".format(
                            ev_row[0], match_coverage, match_cov_cutoff))
                        continue

                if prepend_text is None:
                    annot.product_name = blast_annot.product_name
                else:
                    annot.product_name = "{0} {1}".format(prepend_text, blast_annot.product_name)

                if append_text is not None:
                        annot.product_name = "{0} {1}".format(annot.product_name, append_text)

                log_fh.write("INFO: {1}: Set product name to '{0}' from {3} hit to {2}\n".format(
                        annot.product_name, id, ev_row[0], label))
                
                annot.gene_symbol = blast_annot.gene_symbol
                log_fh.write("INFO: {1}: Set gene_symbol to '{0}' from {3} hit to {2}\n".format(
                        annot.gene_symbol, id, ev_row[0], label))

                for go_annot in blast_annot.go_annotations:
                    annot.add_go_annotation(go_annot)

                for ec_num in blast_annot.ec_numbers:
                    annot.add_ec_number(ec_num)

                # If we get this far we've assigned annotation and don't want to look at any more
                break

    ev_curs.close()
    index_acc_curs.close()
    index_ec_curs.close()
    index_go_curs.close()
                
    
def apply_hmm_evidence(polypeptides=None, ev_conn=None, config=None, ev_config=None, label=None, index_conn=None, log_fh=None):
    """
    Uses HMM evidence to assign functional evidence to polypeptides.  Description of arguments:
    polypeptides: a dict of biothings.Polypeptides objects, keyed on ID
    ev_conn: SQLite3 connection to the parsed HMM search evidence db for this set of searches
    config: The yaml object for the parsed annotation config file
    ev_config: The parsed evidence section for this label within the annotation config file
    label: Label for the evidence track entry within thet annotation config file
    index_conn:  SQLite3 connection to the reference index for the database searched
    """
    default_product = config['general']['default_product_name']
    ev_curs = ev_conn.cursor()
    ev_qry = "SELECT hmm_accession, total_score FROM hmm_hit WHERE qry_id = ? ORDER BY total_hit_eval ASC"

    acc_main_curs = index_conn.cursor()
    acc_main_curs.execute("begin")
    hmm_class_limit = None

    go_curs = index_conn.cursor()
    go_qry = "SELECT go_id FROM hmm_go WHERE hmm_id = ?"

    ec_curs = index_conn.cursor()
    ec_qry = "SELECT ec_id FROM hmm_ec WHERE hmm_id = ?"
    
    if 'class' in ev_config:
        hmm_class_limit = ev_config['class']

    prepend_text = None
    if 'prepend_text' in ev_config:
        prepend_text = ev_config['prepend_text']

    append_text = None
    if 'append_text' in ev_config:
        append_text = ev_config['append_text']

    print("DEBUG: Applying HMM results to {0} polypeptides".format(len(polypeptides)))

    if 'debugging_polypeptide_limit' in config['general']:
        DEBUG_LIMIT = config['general']['debugging_polypeptide_limit']
    
    for id in polypeptides:
        #print("DEBUG: Parsing {0} evidence for polypeptide ID {1}".format(label, id))
        polypeptide = polypeptides[id]
        annot = polypeptide.annotation

        DEBUG_LIMIT = DEBUG_LIMIT - 1
        if DEBUG_LIMIT == 0:
            break

        if config['general']['allow_attributes_from_multiple_sources'] == 'Yes':
            raise Exception("ERROR: Support for the general:allow_attributes_from_multiple_sources=Yes setting not yet implemented")
        else:
            # If this has changed already, it has already been annotated
            if annot.product_name != default_product: continue

            for ev_row in ev_curs.execute(ev_qry, (polypeptide.id,)):
                acc_main_qry = "SELECT version, hmm_com_name, ec_num, isotype, id FROM hmm WHERE (version = ? or accession = ?)"

                if hmm_class_limit is None:
                    acc_main_qry_args = (ev_row[0], ev_row[0], )
                else:
                    acc_main_qry += " AND isotype = ?"
                    acc_main_qry_args = (ev_row[0], ev_row[0], hmm_class_limit)
                
                for acc_main_row in acc_main_curs.execute(acc_main_qry, acc_main_qry_args):
                    if prepend_text is None:
                        annot.product_name = acc_main_row[1]
                    else:
                        annot.product_name = "{0} {1}".format(prepend_text, acc_main_row[1])

                    if append_text is not None:
                        annot.product_name = "{0} {1}".format(annot.product_name, append_text)
                        
                    log_fh.write("INFO: {1}: Set product name to '{0}' from {3} hit to {2}, isotype:{3}\n".format(
                        annot.product_name, id, ev_row[0], hmm_class_limit, label))

                    if acc_main_row[2] is not None:
                        annot.gene_symbol  = acc_main_row[2]
                        log_fh.write("INFO: {1}: Set gene_symbol to '{0}' from {3} hit to {2}, isotype:{3}\n".format(
                            annot.gene_symbol, id, ev_row[0], hmm_class_limit, label))

                    ## add any matching GO terms
                    for go_row in go_curs.execute(go_qry, (acc_main_row[4],)):
                        annot.add_go_annotation(biocode.annotation.GOAnnotation(go_id=go_row[0]))

                    ## add any matching EC numbers
                    for ec_row in ec_curs.execute(ec_qry, (acc_main_row[4],)):
                        annot.add_ec_number(biocode.annotation.ECAnnotation(number=ec_row[0]))
                    
                break
        
    acc_main_curs.close()
    ev_curs.close()
    go_curs.close()
    ec_curs.close()
        

def apply_lipoprotein_motif_evidence(polypeptides=None, ev_conn=None, config=None, ev_config=None, label=None, log_fh=None):
    default_product = config['general']['default_product_name']
    lipoprotein_motif_default_product = ev_config['product_name']

    ev_curs = ev_conn.cursor()
    ev_qry = """
       SELECT hit_acc, hit_desc, start, stop
         FROM lipoprotein_motif_hit
        WHERE qry_id = ?
    """
    if 'debugging_polypeptide_limit' in config['general']:
        DEBUG_LIMIT = config['general']['debugging_polypeptide_limit']

    print("DEBUG: Applying lipoprotein_motif results to {0} polypeptides".format(len(polypeptides)))

    for id in polypeptides:
        #print("DEBUG: Parsing {0} evidence for polypeptide ID {1}".format(label, id))
        polypeptide = polypeptides[id]
        annot = polypeptide.annotation

        DEBUG_LIMIT = DEBUG_LIMIT - 1
        if DEBUG_LIMIT == 0:
            break

        if config['general']['allow_attributes_from_multiple_sources'] == 'Yes':
            raise Exception("ERROR: Support for the general:allow_attributes_from_multiple_sources=Yes setting not yet implemented")
        else:
            if annot.product_name != default_product: continue

            for ev_row in ev_curs.execute(ev_qry, (polypeptide.id,)):
                log_fh.write("INFO: {0}: Set product name to '{1}' because it had a lipoprotein_motif match to accession:{2}, description:{3}\n".format(
                            id, lipoprotein_motif_default_product, ev_row[0], ev_row[1]))
                annot.product_name = lipoprotein_motif_default_product
                break

    ev_curs.close()


def apply_tmhmm_evidence(polypeptides=None, ev_conn=None, config=None, ev_config=None, label=None, log_fh=None):
    default_product = config['general']['default_product_name']
    tmhmm_default_product = ev_config['product_name']
    min_helical_spans = int(ev_config['min_helical_spans'])
    
    ev_curs = ev_conn.cursor()
    ev_qry = """
       SELECT th.id, count(tp.hit_id)
         FROM tmhmm_hit th
              JOIN tmhmm_path tp ON th.id=tp.hit_id
        WHERE th.qry_id = ?
    """
    if 'debugging_polypeptide_limit' in config['general']:
        DEBUG_LIMIT = config['general']['debugging_polypeptide_limit']

    print("DEBUG: Applying TMHMM results to {0} polypeptides".format(len(polypeptides)))

    for id in polypeptides:
        #print("DEBUG: Parsing {0} evidence for polypeptide ID {1}".format(label, id))
        polypeptide = polypeptides[id]
        annot = polypeptide.annotation

        DEBUG_LIMIT = DEBUG_LIMIT - 1
        if DEBUG_LIMIT == 0:
            break

        if config['general']['allow_attributes_from_multiple_sources'] == 'Yes':
            raise Exception("ERROR: Support for the general:allow_attributes_from_multiple_sources=Yes setting not yet implemented")
        else:
            if annot.product_name != default_product: continue

            for ev_row in ev_curs.execute(ev_qry, (polypeptide.id,)):
                if ev_row[1] >= min_helical_spans:
                    annot.product_name = tmhmm_default_product
                    log_fh.write("INFO: {0}: Set product name to '{1}' because it had a TMHMM prediction of {2} transmembrane helices\n".format(
                            id, tmhmm_default_product, ev_row[1]))
                    annot.add_go_annotation( biocode.annotation.GOAnnotation(go_id='0016021') )
                    break

    ev_curs.close()


def check_configuration(conf, userargs):
    """
    Performs any basic checks on the annotation configuration file format/syntax/values.  Ideally done
    before most of the rest of the script to save wasted compute time.
    """
    # make sure each of the expected sections are there
    for section in ['general', 'indexes', 'input', 'order', 'evidence']:
        if section not in conf:
            raise Exception("ERROR: Expected a section called '{0}' in the annotation config file, but didn't find one.".format(section))

    # make sure the input section has at least fasta defined
    if 'polypeptide_fasta' not in conf['input']:
        raise Exception("ERROR: You must at least define 'polypeptide_fasta' data in the 'input' section of the annotation config file")

    # the output format should be one of the recognized ones.
    supported_output_formats = ['fasta', 'gff3']
    userargs.output_format = userargs.output_format.lower()
    #userargs_output_format = userargs.output_format.split('/')
    #userargs_file = userargs_output_format[-1]
    #userargs_file = userargs_file.split('.')
    #if userargs_file[-1] not in supported_output_formats:
    if userargs.output_format not in supported_output_formats:
        raise Exception("ERROR: The output format specified '{0}' isn't supported.  Please choose from: {1}".format(userargs.output_format,
                                                                                                                    ", ".join(supported_output_formats)))
    # user must have defined input GFF if they've requested GFF as output.
    if userargs.output_format == 'gff3':
        if 'gff3' not in conf['input'] or conf['input']['gff3'] is None or len(conf['input']['gff3']) == 0:
            raise Exception("ERROR:  If you requested gff3 formatted output, you must specify the corresponding gff3 input file (which provides the genomic coordinates of the features involved.)")
    
    # make sure there aren't any indexes referenced in the evidence section which are not defined in the indexes section
    indexes = list()
    for label in conf['indexes']:
        indexes.append(label)

    for item in conf['evidence']:
        if 'index' in item and item['index'] not in indexes:
            raise Exception("ERROR: Evidence item '{0}' references an index '{1}' not found in the indexes section of the config file".format(item['label'], item['index']))

    # Any indexes defined should actually exist
    for label in conf['indexes']:
        if not os.path.exists(conf['indexes'][label]):
            raise Exception("ERROR: Index with label:{0} and path:{1} couldn't be found".format(label, conf['indexes'][label]))

def get_sprot_accession_info(conn=None, accession=None, config=None):
    curs = conn.cursor()
    ec_curs = conn.cursor()
    go_curs = conn.cursor()
    annot = biocode.annotation.FunctionalAnnotation(product_name=config['general']['default_product_name'])

    # First we need to get the ID from the accession
    qry = """
          SELECT u.id, ua.accession, u.full_name, u.symbol, ua.res_length
          FROM uniprot_sprot u
               JOIN uniprot_sprot_acc ua ON u.id=ua.id
          WHERE ua.accession = ?
          """

    # accession is in the format: sp|A4YVG3|RRF_BRASO
    abbrev, sprot_acc, sprot_id = accession.split('|')
    
    for row in curs.execute(qry, (sprot_acc,)):
        sprot_id = row[0]
        annot.product_name = row[2]
        annot.gene_symbol = row[3]
        annot.other_attributes['ref_len'] = row[4]

        qry = "SELECT ec_num FROM uniprot_sprot_ec WHERE id = ?"
        for ec_row in ec_curs.execute(qry, (sprot_id,)):
            annot.add_ec_number( biocode.annotation.ECAnnotation(number=ec_row[0]) )

        qry = "SELECT go_id FROM uniprot_sprot_go WHERE id = ?"
        for go_row in go_curs.execute(qry, (sprot_id,)):
            annot.add_go_annotation( biocode.annotation.GOAnnotation(go_id=go_row[0], with_from=row[1]) )

    curs.close()
    ec_curs.close()
    go_curs.close()
    return annot

def get_uniref_accession_info(conn=None, accession=None, config=None, acc_curs=None, ec_curs=None, go_curs=None):
    annot = biocode.annotation.FunctionalAnnotation(product_name=config['general']['default_product_name'])

    if accession.startswith('UniRef100_'):
        accession = accession.lstrip('UniRef100_')

    # First we need to get the ID from the accession
    qry = """
          SELECT u.id, ua.accession, u.full_name, u.symbol, ua.res_length
          FROM uniref u
               JOIN uniref_acc ua ON u.id=ua.id
          WHERE ua.accession = ?
          """

    for row in acc_curs.execute(qry, (accession,)):
        uniref_id = row[0]
        annot.product_name = row[2]
        annot.gene_symbol = row[3]
        annot.other_attributes['ref_len'] = row[4]
        annot.other_attributes['is_characterized'] = row[5]

        qry = "SELECT ec_num FROM uniref_ec WHERE id = ?"
        for ec_row in ec_curs.execute(qry, (uniref_id,)):
            annot.add_ec_number( biocode.annotation.ECAnnotation(number=ec_row[0]) )

        qry = "SELECT go_id FROM uniref_go WHERE id = ?"
        for go_row in go_curs.execute(qry, (uniref_id,)):
            annot.add_go_annotation( biocode.annotation.GOAnnotation(go_id=go_row[0], with_from=row[1]) )

    return annot
        
def get_or_create_db_connections(type_ev=None, configuration=None, evidence=None, label=None,
                                 db_conn=None, output_base=None):
    """
    type_ev must be either 'hmm_ev', 'blast_ev' or 'tmhmm_ev'
    """
    if type_ev in ['hmm_ev', 'blast_ev']:
        index_label = evidence[label]['index']
    elif type_ev in ['tmhmm_ev', 'lipoprotein_motif_ev']:
        index_label = None

    # Use any existing index connection, else attach to it.
    index_conn = None
    if index_label is not None:
        if index_label in db_conn:
            index_conn = db_conn[index_label]
        else:
            try:
                index_conn = sqlite3.connect(configuration['indexes'][index_label])
            except sqlite3.OperationalError as e:
                raise Exception("ERROR: Failed to connect to evidence database {0} because {1}".format(configuration['indexes'][index_label], e))
                
            db_conn[index_label] = index_conn

    # Attach to or create an evidence database
    ev_db_path = "{0}.{1}.sqlite3".format(output_base, type_ev)
    if os.path.exists(ev_db_path):
        ev_db_conn = sqlite3.connect(ev_db_path)
    else:
        ev_db_conn = sqlite3.connect(ev_db_path)
        if type_ev == 'hmm_ev':
            initialize_hmm_results_db(ev_db_conn)
        elif type_ev == 'blast_ev':
            initialize_blast_results_db(ev_db_conn)
        elif type_ev == 'tmhmm_ev':
            initialize_tmhmm_results_db(ev_db_conn)
        elif type_ev == 'lipoprotein_motif_ev':
            initialize_lipoprotein_motif_results_db(ev_db_conn)

    db_conn[type_ev] = ev_db_conn

    # only parse the evidence if the list isn't already in the database
    if not already_indexed(path=evidence[label]['path'], index=ev_db_conn):

        if type_ev == 'hmm_ev':
            index_hmmer3_htab(path=evidence[label]['path'], index=ev_db_conn)
            ev_db_conn.commit()
            # update the database search indexes
            hmm_ev_db_curs = ev_db_conn.cursor()
            hmm_ev_db_curs.execute("DROP INDEX IF EXISTS hmm_hit__qry_id")
            hmm_ev_db_curs.execute("CREATE INDEX hmm_hit__qry_id ON hmm_hit (qry_id)")
            hmm_ev_db_curs.close()
        elif type_ev == 'blast_ev':
            index_rapsearch2_m8(path=evidence[label]['path'], index=ev_db_conn)
            ev_db_conn.commit()
            blast_ev_db_curs = ev_db_conn.cursor()
            blast_ev_db_curs.execute("DROP INDEX IF EXISTS blast_hit__qry_id")
            blast_ev_db_curs.execute("CREATE INDEX blast_hit__qry_id ON blast_hit (qry_id)")
            blast_ev_db_curs.close()
        elif type_ev == 'tmhmm_ev':
            index_tmhmm_raw(path=evidence[label]['path'], index=ev_db_conn)
            ev_db_conn.commit()
            tmhmm_ev_db_curs = ev_db_conn.cursor()
            tmhmm_ev_db_curs.execute("DROP INDEX IF EXISTS tmhmm_hit__qry_id")
            tmhmm_ev_db_curs.execute("CREATE INDEX tmhmm_hit__qry_id ON tmhmm_hit (qry_id)")
            tmhmm_ev_db_curs.execute("DROP INDEX IF EXISTS tmhmm_path__hit_id")
            tmhmm_ev_db_curs.execute("CREATE INDEX tmhmm_path__hit_id ON tmhmm_path (hit_id)")
            tmhmm_ev_db_curs.close()
        elif type_ev == 'lipoprotein_motif_ev':
            index_lipoprotein_motif(path=evidence[label]['path'], index=ev_db_conn)
            ev_db_conn.commit()
            lipo_ev_db_curs = ev_db_conn.cursor()
            lipo_ev_db_curs.execute("DROP INDEX IF EXISTS lipoprotein_motif_hit__qry_id")
            lipo_ev_db_curs.execute("CREATE INDEX lipoprotein_motif_hit__qry_id ON lipoprotein_motif_hit (qry_id)")
            lipo_ev_db_curs.close()

    ev_db_conn.commit()
    return (index_conn, ev_db_conn)

def get_files_from_path(path):
    # Can pass a single file, list file or comma-separated combination.  We'll rely on the .list file extension here
    path_entries = path.split(',')
    paths = list()
    for path in path_entries:
        if path.endswith('.list'):
            paths.extend(biocode.utils.read_list_file(path))
        else:
            paths.extend([path])

    return paths
        
def index_hmmer3_htab(path=None, index=None):
    curs = index.cursor()
    parsing_errors = 0

    qry = """
          INSERT INTO hmm_hit (qry_id, qry_start, qry_end, hmm_accession, hmm_length, hmm_start, hmm_end, 
                               domain_score, total_score, total_score_tc, total_score_nc, total_score_gc,
                               domain_score_tc, domain_score_nc, domain_score_gc, total_hit_eval,
                               domain_hit_eval)
          VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
    """

    paths = get_files_from_path(path)

    for file in paths:
        print("INFO: parsing: {0}".format(file))
        for line in open(file):
            line = line.rstrip()
            cols = line.split("\t")

            ## not sure what this is, but some lines have columns 7+ as these values:
            #   hmm       to      ali     to              bias
            if cols[6] == 'hmm': continue

            # the following columns in the htab files are nullable, which are filled with various widths of '-' characters
            for i in (6,7,8,9,11,19,20):
                if '--' in cols[i]:
                    cols[i] = None
                elif len(cols[i]) > 0:
                    try:
                        cols[i] = float(cols[i])
                    except ValueError:
                        #print("DEBUG: Could not convert this to a float: {0}".format(cols[i]))
                        parsing_errors += 1
                        cols[i] = None
                
            curs.execute(qry, (cols[5], cols[8], cols[9], cols[0], int(cols[2]), cols[6], cols[7], 
                               cols[11], float(cols[12]), float(cols[17]), float(cols[18]), float(cols[23]),
                               float(cols[21]), float(cols[22]), float(cols[24]), cols[19], cols[20]))

    curs.execute("INSERT INTO data_sources (source_path) VALUES (?)", (path,))
    curs.close()

    if parsing_errors > 0:
        print("WARN: There were {0} parsing errors (columns converted to None) when processing {1}\n".format(parsing_errors, path))

def index_lipoprotein_motif(path=None, index=None):
    curs = index.cursor()
    parsing_errors = 0

    qry = """
          INSERT INTO lipoprotein_motif_hit (qry_id, hit_acc, hit_desc, start, stop)
          VALUES (?, ?, ?, ?, ?)
    """

    paths = get_files_from_path(path)

    # http://www.diveintopython3.net/xml.html
    for file in paths:
        print("INFO: parsing: {0}".format(file))
        tree = etree.parse(file)
        for elem in tree.iterfind('Definitions/Sequences/Sequence'):
            qry_id = elem.attrib['id']

            for feature in elem.iterfind('Feature-tables/Feature-table/Feature'):
                title = feature.attrib['title']
                m = re.match("(\S+) \:\: (.+)", title)
                if m:
                    hit_acc = m.group(1)
                    hit_desc = m.group(2)
                    interval = feature.find('Interval-loc')
                    curs.execute(qry, (qry_id, hit_acc, hit_desc, interval.attrib['startpos'],
                                       interval.attrib['endpos']))
                else:
                    parsing_errors += 1
                    print("WARN: Unable to parse accession and description from title in file: {0}".format(file))

    curs.execute("INSERT INTO data_sources (source_path) VALUES (?)", (path,))
    curs.close()
                
        
def index_rapsearch2_m8(path=None, index=None):
    curs = index.cursor()
    parsing_errors = 0

    # The E-value column can be either the E-value directly or log(E-value), depending on
    #  the version and options used.  Luckily, the header line tells us which it is.
    logged_eval = False

    qry = """
          INSERT INTO blast_hit (qry_id, sbj_id, align_len, qry_start, qry_end, sbj_start,
                                 sbj_end, perc_identity, eval, bit_score)
          VALUES (?,?,?,?,?,?,?,?,?,?)
    """

    paths = get_files_from_path(path)
    
    for file in paths:
        print("INFO: parsing: {0}".format(file))
        for line in open(file):
            if line[0] == '#':
                m = re.search('log\(e\-value\)', line)
                if m:
                    logged_eval = True
                continue
            
            line = line.rstrip()
            cols = line.split("\t")

            if len(cols) != 12:
                continue

            if logged_eval == True:
                try:
                    cols[10] = math.pow(10, float(cols[10]))
                except OverflowError:
                    # RapSearch2 sometimes reports miniscule e-values, such a log(eval) of > 1000
                    #  These are outside of the range of Python's double.  In my checking of these
                    #  though, their alignments don't warrant a low E-value at all.  Skipping them.
                    print("WARN: Skipping a RAPSearch2 row:  overflow error converting E-value ({0}) on line: {1}".format(cols[10], line))
                    parsing_errors += 1
                    continue
                
            curs.execute(qry, (cols[0], cols[1], int(cols[3]), int(cols[6]), int(cols[7]), int(cols[8]),
                               int(cols[9]), float(cols[2]), cols[10], float(cols[11])))

    curs.execute("INSERT INTO data_sources (source_path) VALUES (?)", (path,))
    curs.close()
    if parsing_errors > 0:
        print("WARN: There were {0} parsing errors (match rows skipped) when processing {1}\n".format(parsing_errors, path))


def index_tmhmm_raw(path=None, index=None):
    """
    Notes from the esteemed M Giglio:
    The GO term to use would be GO:0016021 "integral component of membrane"
    Or if you want to be more conservative you could go with GO:0016020 "membrane"
    
    Depends on the evidence. For the prok pipe we are pretty conservative, we require five TMHMM
    domains and then we call it putative integral membrane protein. 
    On ECO - in fact Marcus and I are the developers of ECO.  It is an ontology of evidence types.
    An annotation to an ECO term is used in conjunction with another annotation, like a GO term
    (but many other types of annotation can, and are, used with ECO). It provides additional
    information about the annotation. In fact for GO, the assignment of an evidence term along
    with a GO term is a required part of a GO annotation. (ECO terms are the "evidence codes" in GO.)
    INPUT: Expected TMHMM input (all HTML lines are skipped)
    # CHARM010_V2.mRNA.887 Length: 904
    # CHARM010_V2.mRNA.887 Number of predicted TMHs:  6
    # CHARM010_V2.mRNA.887 Exp number of AAs in TMHs: 133.07638
    # CHARM010_V2.mRNA.887 Exp number, first 60 AAs:  21.83212
    # CHARM010_V2.mRNA.887 Total prob of N-in:        0.99994
    # CHARM010_V2.mRNA.887 POSSIBLE N-term signal sequence
    CHARM010_V2.mRNA.887	TMHMM2.0	inside	     1    11
    CHARM010_V2.mRNA.887	TMHMM2.0	TMhelix	    12    34
    CHARM010_V2.mRNA.887	TMHMM2.0	outside	    35   712
    CHARM010_V2.mRNA.887	TMHMM2.0	TMhelix	   713   735
    CHARM010_V2.mRNA.887	TMHMM2.0	inside	   736   755
    CHARM010_V2.mRNA.887	TMHMM2.0	TMhelix	   756   773
    CHARM010_V2.mRNA.887	TMHMM2.0	outside	   774   782
    CHARM010_V2.mRNA.887	TMHMM2.0	TMhelix	   783   805
    CHARM010_V2.mRNA.887	TMHMM2.0	inside	   806   809
    CHARM010_V2.mRNA.887	TMHMM2.0	TMhelix	   810   832
    CHARM010_V2.mRNA.887	TMHMM2.0	outside	   833   871
    CHARM010_V2.mRNA.887	TMHMM2.0	TMhelix	   872   894
    CHARM010_V2.mRNA.887	TMHMM2.0	inside	   895   904
    """
    curs = index.cursor()
    parsing_errors = 0

    hit_qry = """
       INSERT INTO tmhmm_hit (qry_id, tmh_count, num_aa_in_tmhs, num_aa_in_f60, total_prob_n_in)
       VALUES (?, ?, ?, ?, ?)
    """

    path_qry = """
       INSERT INTO tmhmm_path (hit_id, locus, start, stop)
       VALUES (?, ?, ?, ?)
    """

    paths = get_files_from_path(path)

    for file in paths:
        print("INFO: parsing: {0}".format(file))
        last_qry_id = None
        current_hit_id = None
        current_path = list()
        tmh_count = num_aa_in_tmhs = num_aa_in_f60 = total_prob_n_in = 0
        
        for line in open(file):
            # skip the HTML lines
            if line.startswith('<'): continue

            m = Match()
            if m.match("# (.+?)\s+Length: \d+", line): # this line marks a new result
                current_id = m.m.group(1)

                # purge the previous result
                if last_qry_id is not None:
                    curs.execute(hit_qry, (last_qry_id, tmh_count, num_aa_in_tmhs, num_aa_in_f60, total_prob_n_in))
                    current_hit_id = curs.lastrowid

                for span in current_path:
                    curs.execute(path_qry, (current_hit_id, span[2], int(span[3]), int(span[4])))

                # reset
                last_qry_id = current_id
                current_helix_count = tmh_count = num_aa_in_tmhs = num_aa_in_f60 = total_prob_n_in = 0
                current_path = list()
                
            elif m.match(".+Number of predicted TMHs:\s+(\d+)", line):
                tmh_count = int(m.m.group(1))
            elif m.match(".+Exp number of AAs in TMHs:\s+([0-9\.]+)", line):
                num_aa_in_tmhs = float(m.m.group(1))
            elif m.match(".+Exp number, first 60 AAs:\s+([0-9\.]+)", line):
                num_aa_in_f60 = float(m.m.group(1))
            elif m.match(".+Total prob of N-in:\s+([0-9\.]+)", line):
                total_prob_n_in = float(m.m.group(1))
            else:
                if line[0] == '#': continue
                
                cols = line.split()
                if len(cols) == 5:
                    current_path.append(cols)

        # don't forget to do the last entry
        curs.execute(hit_qry, (last_qry_id, tmh_count, num_aa_in_tmhs, num_aa_in_f60, total_prob_n_in))
        current_hit_id = curs.lastrowid
        for span in current_path:
            curs.execute(path_qry, (current_hit_id, span[2], int(span[3]), int(span[4])))
        

    curs.execute("INSERT INTO data_sources (source_path) VALUES (?)", (path,))
    curs.close()
    index.commit()
    
    if parsing_errors > 0:
        print("WARN: There were {0} parsing errors (match rows skipped) when processing {1}\n".format(parsing_errors, path))

    
def initialize_blast_results_db(conn):
    curs = conn.cursor()

    curs.execute("""
        CREATE TABLE blast_hit (
            id                integer primary key,
            qry_id            text,
            sbj_id            text,
            align_len         integer,
            qry_start         integer,
            qry_end           integer,
            sbj_start         integer,
            sbj_end           integer,
            perc_identity     real,
            eval              real,
            bit_score         real
        )
    """)

    curs.execute("""
        CREATE TABLE data_sources (
            id                integer primary key,
            source_path         text
        )
    """)

    curs.close()
    conn.commit()
                
def initialize_hmm_results_db(conn):
    curs = conn.cursor()

    curs.execute("""
        CREATE TABLE hmm_hit (
            id                integer primary key,
            qry_id            text,
            qry_start         integer,
            qry_end           integer,
            hmm_accession     text,
            hmm_length        integer,
            hmm_start         integer,
            hmm_end           integer,
            domain_score      real,
            total_score       real,
            total_score_tc    real,
            total_score_nc    real,
            total_score_gc    real,
            total_hit_eval    real,
            domain_score_tc   real,
            domain_score_nc   real,
            domain_score_gc   real,
            domain_hit_eval   real
        )
    """)

    curs.execute("""
        CREATE TABLE data_sources (
            id                integer primary key,
            source_path         text
        )
    """)

    curs.close()
    conn.commit()

def initialize_lipoprotein_motif_results_db(conn):
    curs = conn.cursor()

    curs.execute("""
        CREATE TABLE lipoprotein_motif_hit (
            id                integer primary key,
            qry_id            text,
            hit_acc           text,
            hit_desc          text,
            start             integer,
            stop              integer
        )
    """)

    curs.execute("""
        CREATE TABLE data_sources (
            id                integer primary key,
            source_path         text
        )
    """)

    curs.close()
    conn.commit()

def initialize_tmhmm_results_db(conn):
    curs = conn.cursor()

    """
    tmh_count: Number of predicted Trans-Membrane Helices (TMHs)
    num_aa_in_tmhs: Exp number of AAs in TMHs
    num_aa_in_f60: Exp number of AAs in TMHs within first 60 AAs
    
    """
    curs.execute("""
        CREATE TABLE tmhmm_hit (
            id                integer primary key,
            qry_id            text,
            tmh_count         float,
            num_aa_in_tmhs    float,
            num_aa_in_f60     float,
            total_prob_n_in   float
        )
    """)

    curs.execute("""
        CREATE TABLE tmhmm_path (
            hit_id            integer,
            locus             text,
            start             integer,
            stop              integer
        )
    """)

    curs.execute("""
        CREATE TABLE data_sources (
            id                integer primary key,
            source_path         text
        )
    """)

    curs.close()
    conn.commit()

        
def initialize_polypeptides( log_fh, fasta_file, default_name ):
    '''
    Reads a FASTA file of (presumably) polypeptide sequences and creates a dict of Polypeptide
    objects, keyed by ID, with bioannotation.FunctionalAnnotation objects attached.
    '''
    seqs = biocode.utils.fasta_dict_from_file( fasta_file )

    polypeptides = dict()

    for seq_id in seqs:
        polypeptide = biocode.things.Polypeptide( id=seq_id, length=len(seqs[seq_id]['s']), residues=seqs[seq_id]['s'] )
        annotation = biocode.annotation.FunctionalAnnotation(product_name=default_name)
        log_fh.write("INFO: {0}: Set initial product name to '{1}'\n".format(seq_id, default_name))
        polypeptide.annotation = annotation
        
        polypeptides[seq_id] = polypeptide
    
    return polypeptides
        
def parse_evidence_config(conf):
    """
    Parses the 'evidence' section of the annotation config file, and returns a dict where each key
    is the label of that evidence and the value is a dict of the other key/value pairs
    """
    ev = dict()

    for entry in conf['evidence']:
        # make sure there aren't duplicates
        label = entry['label']

        if label in ev:
            raise Exception("ERROR: duplicate label found in evidence track: {0}".format(label))
        else:
            ev[label] = dict()
            for key in entry:
                if key is not 'label':
                    ev[label][key] = entry[key]

    return ev


def perform_final_checks(polypeptides=None, config=None, log_fh=None):
    """
    Does a round of checks we want to perform on an annotated set of polypeptides
    before exporting them.  Currently:
    - Make sure a gene product name is assigned.  This might accidentally become "None" if
      a match existed to a subject which didn't have a name properly entered in the index.
    """
    for id in polypeptides:
        polypeptide = polypeptides[id]

        if polypeptide.annotation.product_name is None:
            log_fh.write("WARNING: {0}: Somehow made it through annotation with no product name.  Setting to default\n".format(id))
            polypeptide.annotation.product_name = config['general']['default_product_name']
    
    
class Match(object):
    """
    Python doesn't really have a good syntax for doing a series of conditional checks on a line if you want
    to also capture part of the matter and use it.  This wraps the match object so that we can do exactly that.
    Example use:
    match = Match() 
    if match.match(pattern1, string1): 
       do_something( print(match.m.group(1)) ) 
    elif match.match(pattern2, string2): 
       do_something_else( print(match.m.group(1)) ) 
    """
    def __init__(self): 
        self.m = None 

    def match(self, *args, **kwds): 
        self.m = re.match(*args, **kwds)
        return self.m is not None 

if __name__ == '__main__':
    cProfile.run('main()', "{0}/FALCON.profile.{1}.dat".format(out_dir, os.getpid()))
