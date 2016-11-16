#!/usr/bin/env python
import argparse
import gzip
import collections
import os
import logging

logger = logging.getLogger(os.path.basename(__file__))
logging.basicConfig(level=logging.INFO)

__author__ = 'Alexander Junge (alexander.junge@gmail.com), RTH, University of Copenhagen'


def extant_file(x):
    """
    'Type' for argument parsing - checks that file exists but does not open.
    :param x: file to be tested for existence
    """
    if not os.path.exists(x):
        raise argparse.ArgumentError(x, "{0} does not exist".format(x))
    return x


parser = argparse.ArgumentParser(description="""
                                 Creates RNA alias file to be placed on website and used in identifier resolution.
                                 File is similar to STRING's protein.aliases.txt.gz file:
                                 http://string-db.org/newstring_download/protein.aliases.v10.txt.gz""")
parser.add_argument('aggregated_master_file', help='file containing interactions collected from all master files',
                    type=extant_file)
parser.add_argument('tax_identifiers', type=str, nargs='+', help='taxonomy identifiers of organisms to be used')
parser.add_argument('--print_universe', action='store_true', help='if set, aggregated_master_file & tax_identifiers '
                                                                  'arguments are '
                                                                  'ignored and all known aliases are returned.')
args = parser.parse_args()

tax_identifiers = args.tax_identifiers
MASTER_FILE_PATH = args.aggregated_master_file
print_universe = args.print_universe

logger.info("Creating aliases file for taxonomy identifiers: {}".format(', '.join(tax_identifiers)))

NCRNA_ALIAS_FILE_PATH = 'id_dictionaries/ncRNAaliasfile.tsv.gz'
MIRNA_ALIAS_FILE_PATH = 'id_dictionaries/miRBase_20_mir_aliases.tsv'

tax_id_to_existing_nodes = collections.defaultdict(set)
with open(MASTER_FILE_PATH, 'r') as master_file:
    for line in master_file:
        tax_id, node1, node2, directed, evidence, score, source, link_out, comment = line.rstrip('\n').split('\t')
        if tax_id in tax_identifiers:
            tax_id_to_existing_nodes[tax_id].add(node1)
            tax_id_to_existing_nodes[tax_id].add(node2)
logger.info('Found %d nodes in aggregated master file %s. Considered %d taxonomy identifiers:\n%s' %
             (sum([len(v) for v in tax_id_to_existing_nodes.itervalues()]), MASTER_FILE_PATH,
              len(tax_identifiers), ' '.join(tax_identifiers)))

rna_tax_identifier_to_tax_alias_sources = collections.defaultdict(list)
with open(MIRNA_ALIAS_FILE_PATH, 'r') as mirna_alias_file:
    for line in mirna_alias_file:
        tax_id, identifier, alias, priority = line.rstrip('\n').split('\t')
        source = 'miRBase'
        tax_alias_source = [int(tax_id), alias, {source}]
        key = (int(tax_id), identifier)

        if not print_universe and identifier not in tax_id_to_existing_nodes[tax_id]:
            continue

        if key in rna_tax_identifier_to_tax_alias_sources:
            alias_added = False
            for prev_tax_alias_sources in rna_tax_identifier_to_tax_alias_sources[key]:
                if prev_tax_alias_sources[0] == int(tax_id) and prev_tax_alias_sources[1] == alias:
                    prev_tax_alias_sources[2].add(source)
                    alias_added = True
            if not alias_added:
                rna_tax_identifier_to_tax_alias_sources[key].append(tax_alias_source)
        else:
            rna_tax_identifier_to_tax_alias_sources[key].append(tax_alias_source)
logger.info('Extracted %d identifiers and %d aliases from alias file %s.' %
             (len(rna_tax_identifier_to_tax_alias_sources),
              sum([len(v) for v in rna_tax_identifier_to_tax_alias_sources.itervalues()]),
              MIRNA_ALIAS_FILE_PATH))

with gzip.open(NCRNA_ALIAS_FILE_PATH, 'r') as ncrna_alias_file:
    for line in ncrna_alias_file:
        tax_id, identifier, alias, source = line.rstrip('\n').split('\t')

        tax_alias_source = [int(tax_id), alias, {source}]
        key = (int(tax_id), identifier)

        if not print_universe and identifier not in tax_id_to_existing_nodes[tax_id]:
            continue

        if key in rna_tax_identifier_to_tax_alias_sources:
            alias_added = False
            for prev_tax_alias_sources in rna_tax_identifier_to_tax_alias_sources[key]:
                if prev_tax_alias_sources[0] == int(tax_id) and prev_tax_alias_sources[1] == alias:
                    prev_tax_alias_sources[2].add(source)
                    alias_added = True
            if not alias_added:
                rna_tax_identifier_to_tax_alias_sources[key].append(tax_alias_source)
        else:
            rna_tax_identifier_to_tax_alias_sources[key].append(tax_alias_source)
logger.info('Extracted %d identifiers and %d aliases from alias files %s and %s.' %
             (len(rna_tax_identifier_to_tax_alias_sources),
              sum([len(v) for v in rna_tax_identifier_to_tax_alias_sources.itervalues()]),
              MIRNA_ALIAS_FILE_PATH,
              NCRNA_ALIAS_FILE_PATH))

# print header
header_row = '## species_ncbi_taxon_id ## ncrna_id ## alias ## source ##'
print(header_row)

# print sorted by taxonomy id and RNA identifier
sorted_rna_tax_id_keys = list(rna_tax_identifier_to_tax_alias_sources.keys())
sorted_rna_tax_id_keys.sort()
for rna_tax_id in sorted_rna_tax_id_keys:
    tax_id, rna_id = rna_tax_id
    rna_tax_alias_sources_list = rna_tax_identifier_to_tax_alias_sources[rna_tax_id]
    for rna_tax_alias_sources in rna_tax_alias_sources_list:
        tax, rna_alias, sources = rna_tax_alias_sources
        if tax_id != tax:
            raise IOError("Identifier: " + str(rna_tax_id) + " , Alias: " + str((tax, rna_alias)))
        sources_string = ''
        for i, src in enumerate(sources):
            src_underscores = src.replace(' ', '_')
            if i == 0:
                delimit = ''
            else:
                delimit = ' '
            sources_string += delimit + src_underscores
        row = (str(tax), rna_id, rna_alias, sources_string)
        print('\t'.join(row))
logger.info('Done.' + os.linesep + os.linesep)