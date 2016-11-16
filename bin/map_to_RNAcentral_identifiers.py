#!/usr/bin/env python
import argparse
import gzip
import logging
import os
import urllib
from collections import defaultdict

import re

from bin import stringrnautils

logging.basicConfig(level=logging.INFO)

__author__ = 'Alexander Junge (alexander.junge@gmail.com)'


def is_existing_dir(arg):
    if not os.path.isdir(arg):
        parser.error('Given path %s does not point to an existing directory. Exiting!')
    return arg


parser = argparse.ArgumentParser(description="""
                                 Maps miRNA and other ncRNA identifiers used in RAIN to RNAcentral identifier.
                                 """)
parser.add_argument('-r', '--rnacentral_version', default='4.0', help='RNAcentral version to be used for mapping')
parser.add_argument('-d', '--data_dir', default='data', type=is_existing_dir,
                    help='directory the RNAcentral identifier mapping file is downloaded to, directory must exist '
                         'prior to execution')
parser.add_argument('-o', '--output_mapping_file', default='RNAcentral_mapping.tsv.gz',
                    help='The file name of the RNAcentral mapping file created by this script.')
parser.add_argument('-i', '--id_dictionaries_directory', default='id_dictionaries',
                    help='The directory the RNAcentral mapping file will be written to.')

args = parser.parse_args()
rnacentral_version = args.rnacentral_version
# rnacentral_version = '4.0'
data_dir = args.data_dir
# data_dir = 'data'
output_path = os.path.join(args.id_dictionaries_directory, args.output_mapping_file)
# output_path = os.path.join('id_dictionaries', 'RNAcentral_mapping.tsv.gz')

rnacentral_uri = 'ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/%s/id_mapping/id_mapping.tsv.gz' % \
                 rnacentral_version
rnacentral_file = os.path.join(data_dir, 'RNAcentral_%s_id_mapping.tsv.gz' % rnacentral_version)
if not os.path.exists(rnacentral_file):
    urllib.urlretrieve(rnacentral_uri, rnacentral_file)

logging.info('Loading ncRNA alias mapper.')
tax_to_src_to_ncrna_alias_to_id = stringrnautils.get_non_coding_rna_alias_mapper_including_source()

logging.info('Loading miRNA alias mapper.')
mirna_alias_to_id = stringrnautils.get_unique_mir_mapper()

logging.info('Loading miRNA taxonomy mapper.')
mirna_id_to_tax = stringrnautils.get_mir_id_to_tax_id_mapper()

id_source_databases = [src for k in tax_to_src_to_ncrna_alias_to_id
                       for src in tax_to_src_to_ncrna_alias_to_id[k].keys()]
id_source_databases.append('miRBase')
id_source_databases = set([src.upper() for src in id_source_databases])
logging.info('ncRNA/miRNA aliases came from %d different sources.' % len(id_source_databases))

tax_to_src_to_id_to_rnacentral_ids = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
line_number = 0
id_sources_not_observed = id_source_databases.copy()
logging.info('Loading RNAcentral identifier mappings from file %s.' % rnacentral_file)
with gzip.open(rnacentral_file) as file_handle:
    # Example line
    # URS0000000001   ENA     GU792786.1:1..200:rRNA  77133
    for line in file_handle:
        line_number += 1
        if line_number % 1000000 == 0:
            logging.debug('Read %d lines in RNAcentral mapping file.' % line_number)
        rnacentral_id, source, source_id, taxonomy_id = line.rstrip('\r\n').split('\t')
        # Maps taxonomy ID used in SGD to 'vanilla' yeast taxon ID
        if taxonomy_id == '559292':
            taxonomy_id = '4932'
        if source not in id_source_databases:
            continue
        if source in id_sources_not_observed:
            id_sources_not_observed.remove(source)
        tax_to_src_to_id_to_rnacentral_ids[taxonomy_id][source][source_id].append(rnacentral_id)
        # if rnacentral_mappings == 100000:  # for testing purposes
        #     break

rnacentral_mapping_count = sum([len(tax_to_src_to_id_to_rnacentral_ids[tax][src][src_id])
                                for tax in tax_to_src_to_id_to_rnacentral_ids
                                for src in tax_to_src_to_id_to_rnacentral_ids[tax]
                                for src_id in tax_to_src_to_id_to_rnacentral_ids[tax][src]])
logging.info('Read %d identifier mappings from file %s.' % (rnacentral_mapping_count, rnacentral_file))

if len(id_sources_not_observed) > 0:
    logging.warn('Of the following %d databases, not a single ID mapping was found in RNAcentral:%s%s' %
                 (len(id_sources_not_observed), os.linesep, ', '.join(id_sources_not_observed)))
    logging.warn('A complete list of ID mapping source data bases is: %s' % ', '.join(id_source_databases))

logging.info('Mapping miRBase aliases to RNAcentral identifiers.')
mirna_aliases_mapped = 0
mirna_aliases_seen = 0
mirna_identifiers_mapped = set()
mirna_identifiers_seen = set()
mirna_id_to_rnacentral_id = defaultdict(set)
for alias, identifier in mirna_alias_to_id.iteritems():
    tax = mirna_id_to_tax[identifier]
    src = 'MIRBASE'
    # Skip aliases that clearly come from precursors
    precursor_accession_re = re.compile('MI\d+')
    precursor_id_re = re.compile('\S+-mir-')
    if precursor_accession_re.match(alias) or precursor_id_re.match(alias):
        continue
    mirna_aliases_seen += 1
    mirna_identifiers_seen.add(identifier)
    if len(tax_to_src_to_id_to_rnacentral_ids[tax][src][alias]) > 0:
        mirna_aliases_mapped += 1
        mirna_identifiers_mapped.add(identifier)
        mirna_id_to_rnacentral_id[identifier].update(tax_to_src_to_id_to_rnacentral_ids[tax][src][alias])
logging.info('%d out of %d miRBase aliases were mapped to RNAcentral identifiers.' % (mirna_aliases_mapped,
                                                                                      mirna_aliases_seen))
logging.info('%d out of %d unique miRBase identifiers were mapped to RNAcentral identifiers.' %
             (len(mirna_identifiers_mapped), len(mirna_identifiers_seen)))
assert len(mirna_identifiers_mapped) == len(mirna_id_to_rnacentral_id)

multi_mapped_mirna_count = 0
for mirna_id, rnacentral_ids in mirna_id_to_rnacentral_id.iteritems():
    if len(rnacentral_ids) > 1:
        multi_mapped_mirna_count += 1
        logging.debug('miRNA %s is mapping to multiple RNAcentral IDs: %s' % (mirna_id, ', '.join(rnacentral_ids)))
if multi_mapped_mirna_count > 0:
    logging.warn('%d out of %d miRBase identifiers were mapped to more than one RNAcentral identifier.' %
                 (multi_mapped_mirna_count, len(mirna_id_to_rnacentral_id)))

logging.info('Mapping remaining ncRNA aliases to RNAcentral identifiers.')
ncrna_aliases_mapped = 0
ncrna_aliases_seen = 0
ncrna_identifiers_mapped = set()
ncrna_identifiers_seen = set()
tax_to_ncrna_id_to_rnacentral_id = defaultdict(lambda: defaultdict(set))
for tax in tax_to_src_to_ncrna_alias_to_id:
    for src in tax_to_src_to_ncrna_alias_to_id[tax]:
        for ncrna_alias, ncrna_id in tax_to_src_to_ncrna_alias_to_id[tax][src].iteritems():
            ncrna_aliases_seen += 1
            ncrna_identifiers_seen.add((tax, ncrna_id))
            if len(tax_to_src_to_id_to_rnacentral_ids[tax][src.upper()][ncrna_alias]) > 0:
                ncrna_aliases_mapped += 1
                ncrna_identifiers_mapped.add((tax, ncrna_id))
                tax_to_ncrna_id_to_rnacentral_id[tax][ncrna_id].update(
                    tax_to_src_to_id_to_rnacentral_ids[tax][src.upper()][ncrna_alias])
logging.info('%d out of %d ncRNA aliases were mapped to RNAcentral identifiers.' % (ncrna_aliases_mapped,
                                                                                    ncrna_aliases_seen))
logging.info('%d out of %d ncRNA identifiers were mapped to RNAcentral identifiers.' % (len(ncrna_identifiers_mapped),
                                                                                        len(ncrna_identifiers_seen)))
multi_mapped_ncrna_count = 0
for tax in tax_to_ncrna_id_to_rnacentral_id:
    for ncrna_id, rnacentral_ids in tax_to_ncrna_id_to_rnacentral_id[tax].iteritems():
        if len(rnacentral_ids) > 1:
            multi_mapped_ncrna_count += 1
        logging.debug('ncRNA %s is mapping to multiple RNAcentral IDs: %s' % (ncrna_id, ', '.join(rnacentral_ids)))
if multi_mapped_ncrna_count > 0:
    logging.warn('%d out of %d ncRNA identifiers were mapped to more than one RNAcentral identifier.' %
                 (multi_mapped_ncrna_count, len(ncrna_identifiers_mapped)))

with gzip.open(output_path, 'w') as output_file:
    for mirna_id, rnacentral_ids in mirna_id_to_rnacentral_id.iteritems():
        for rnacentral_id in rnacentral_ids:
            output_file.write('\t'.join([rnacentral_id, mirna_id_to_tax[mirna_id], mirna_id]))
            output_file.write(os.linesep)
    for tax in tax_to_ncrna_id_to_rnacentral_id:
        for ncrna_id, rnacentral_ids in tax_to_ncrna_id_to_rnacentral_id[tax].iteritems():
            for rnacentral_id in rnacentral_ids:
                output_file.write('\t'.join([rnacentral_id, tax, ncrna_id]))
                output_file.write(os.linesep)
