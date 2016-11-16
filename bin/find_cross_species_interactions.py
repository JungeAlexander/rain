#!/usr/bin/env python
import argparse
import os
import stringrnautils
import logging
logger = logging.getLogger(os.path.basename(__file__))
logging.basicConfig(level=logging.INFO)

__author__ = 'Alexander Junge (alexander.junge@gmail.com), RTH, University of Copenhagen'


def extant_file(x):
    """
    'Type' for argument parsing - checks that file exists but does not open.
    """
    if not os.path.exists(x):
        raise argparse.ArgumentError(x, "{0} does not exist".format(x))
    return x

parser = argparse.ArgumentParser(description=
                                 """Identifies and returns a list of interactions that are suspected to be fawlty
                                    cross-species interactions.
                                    """)
parser.add_argument('rna_aliases_file', help='alias file for all RNAs in STRING-RNA', type=extant_file)
parser.add_argument('master_files', type=extant_file, nargs='+',
                    help='master files to be checked for cross-species interactions')

args = parser.parse_args()

rna_aliases_file = args.rna_aliases_file
master_files = args.master_files

logger.info('Performing sanity check for cross-species interactions.')
string_id_to_aliases_map = stringrnautils.get_string_to_alias_mapper('all', '', '', 10, 'all', True)
rna_id_to_aliases_map = stringrnautils.get_rna_identifiers_in_organism(rna_aliases_file)

not_found_entities = set()
for master_file_path in master_files:
    #logging.info('Reading master file: {}\n'.format(master_file_path))
    with open(master_file_path, 'r') as master_file:
        for line in master_file:
            tax_id, node1, node2, directed, evidence, score, source, link_out, comment = line.rstrip('\n').split('\t')
            org_rna_map = rna_id_to_aliases_map[tax_id]
            org_protein_map = string_id_to_aliases_map[tax_id]
            if not(node1 in org_rna_map or node1 in org_protein_map):
                not_found_entities.add((node1, tax_id))
            if not(node2 in org_rna_map or node2 in org_protein_map):
                not_found_entities.add((node2, tax_id))

if len(not_found_entities) > 0:
    logger.info('{:d} entities were not found in respective protein/ncRNA alias maps.'.format(len(not_found_entities)))
    logger.debug(
        '{:d} entities were not found in respective protein/ncRNA alias maps: {}'.format(len(not_found_entities),
                                                                                         ', '.join(sorted(
                                                                                             [str(x) for x in
                                                                                              not_found_entities]))))

else:
    logger.info('No suspected cross-species interactions found.')
logger.info('Done.' + os.linesep + os.linesep)