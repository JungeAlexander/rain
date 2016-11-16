#!/usr/bin/env python
import argparse
import collections
import stringrnautils
import sys
import os
import logging
logger = logging.getLogger(os.path.basename(__file__))
logging.basicConfig(level=logging.INFO)


__author__ = 'Alexander Junge (alexander.junge@gmail.com)'

parser = argparse.ArgumentParser(description="""
                                 Determines the number of unique interactions per species in STRING-RNA.
                                 Taxonomy IDs of species having the at least a given minimum number of interactions are
                                 printed to stdout.
                                 Prints the number of interactions for all species to stderr.
                                 Furthermore prints a TeX-style table to stderr that lists the number of interactions
                                 of different types (miRNA-mRNA, ncRNA-protein, ncRNA-ncRNA).""")
parser.add_argument('n', type=int,
                    help='min number of interactions required per species to print the species taxonomy ID to stdout.')
parser.add_argument('-id_path',
                    help='path to write id dictionaries',
                    default='id_dictionaries')
parser.add_argument('-in_file',
                    help='path to all.tsv file containing all interactions to be taken into account',
                    default='data/all.tsv')

args = parser.parse_args()
min_interaction_count = args.n
id_dict_dir = args.id_path

logger.info('Filtering species such that only those with a minimum number of {:d} interactions are retained.'.format(
        min_interaction_count))

in_file = args.in_file

tax_to_interactions = collections.defaultdict(set)
with open(in_file, 'r') as file_handle:
    for line in file_handle:
        cols = line.rstrip('\n').split('\t')
        tax, ent1_unsort, ent2_unsort = cols[0:3]
        source = cols[4]
        ent1, ent2 = sorted((ent1_unsort, ent2_unsort))
        tax_to_interactions[tax].add((ent1, ent2, source))

evidence_channel_to_table_name = {'database': 'Curated', 'experimental': 'Experiments', 'prediction': 'Predictions',
                                  'textmining': 'Text mining'
                                  }

all_taxonomy_ids = list(tax_to_interactions.keys())
protein_alias_mapper = stringrnautils.get_string_to_alias_mapper(all_taxonomy_ids, '', '', 10, 'all', True)

mir_alias_mapper = stringrnautils.get_unique_mir_mapper()

# load ncRNA dictionary
ncrna_mapper = stringrnautils.get_non_coding_rna_alias_mapper()  # taxonomy ID -> ncRNA alias -> ncRNA identifier

tax_interaction_count = []
tax_to_miRNA_mRNA_interactions = collections.defaultdict(set)
tax_to_ncRNA_protein_interactions = collections.defaultdict(set)
tax_to_ncRNA_ncRNA_interactions = collections.defaultdict(set)
tax_to_miRNAs = collections.defaultdict(set)
tax_to_ncRNAs = collections.defaultdict(set)
tax_to_proteins = collections.defaultdict(set)
tax_to_source_to_interactions = {}
for tax, interaction_set in tax_to_interactions.iteritems():
    tax_interaction_count.append((tax, len(interaction_set)))
    tax_protein_alias_mapper = protein_alias_mapper[tax]
    tax_miRNAs = tax_to_miRNAs[tax]
    tax_ncRNAs = tax_to_ncRNAs[tax]
    tax_proteins = tax_to_proteins[tax]
    tax_miRNA_mRNA_interactions = tax_to_miRNA_mRNA_interactions[tax]
    tax_ncRNA_protein_interactions = tax_to_ncRNA_protein_interactions[tax]
    tax_ncRNA_ncRNA_interactions = tax_to_ncRNA_ncRNA_interactions[tax]

    logger.info('Processing tax {}.'.format(str(tax)))
    assumed_ncRNA = set()
    ppi_interactions = set()

    for entity1, entity2, source in interaction_set:
        if entity1 < entity2:
            entity_key = (entity1, entity2)
        else:
            entity_key = (entity2, entity1)

        if source not in evidence_channel_to_table_name:
            raise ValueError('Found source {} for interaction, which is not allowed.'.format(source))
        if tax not in tax_to_source_to_interactions:
            tax_to_source_to_interactions[tax] = {}
        if source not in tax_to_source_to_interactions[tax]:
            tax_to_source_to_interactions[tax][source] = set()
        if entity_key in tax_to_source_to_interactions[tax][source]:
            raise ValueError('Interaction of {} and {} occurred multiple times in tax {} evidence channel.'.format(
                entity1, entity2, tax, source))
        tax_to_source_to_interactions[tax][source].add(entity_key)

        ncRNA_mapper_key_entity1 = (tax, entity1)
        entity1_is_protein = entity1 in tax_protein_alias_mapper
        entity1_is_miRNA = entity1 in mir_alias_mapper
        entity1_is_ncRNA = entity1 in ncrna_mapper[tax]
        if entity1_is_ncRNA:
            # used to be protein according to STRING but we consider it a ncRNA
            entity1_is_protein = False

        ncRNA_mapper_key_entity2 = (tax, entity2)
        entity2_is_protein = entity2 in tax_protein_alias_mapper
        entity2_is_miRNA = entity2 in mir_alias_mapper
        entity2_is_ncRNA = entity2 in ncrna_mapper[tax]
        if entity2_is_ncRNA:
            # used to be protein according to STRING but we consider it a ncRNA
            entity2_is_protein = False

        # Count proteins/ncRNAs/miRNAs
        if entity1_is_protein:
            tax_proteins.add(entity1)
        elif entity1_is_ncRNA:
            tax_ncRNAs.add(entity1)
        elif entity1_is_miRNA:
            tax_miRNAs.add(entity1)
        else:
            entity1_is_ncRNA = True
            tax_ncRNAs.add(entity1)
            assumed_ncRNA.add(entity1)

        if entity2_is_protein:
            tax_proteins.add(entity2)
        elif entity2_is_ncRNA:
            tax_ncRNAs.add(entity2)
        elif entity2_is_miRNA:
            tax_miRNAs.add(entity2)
        else:
            entity2_is_ncRNA = True
            tax_ncRNAs.add(entity2)
            assumed_ncRNA.add(entity2)

        if entity1_is_protein:
            if entity2_is_protein:
                ppi_interactions.add(entity1 + '-' + entity2)
                # assuming this is a lncRNA-protein interaction where the lncRNA used to be annotated as a protein
                tax_ncRNA_protein_interactions.add(entity_key)
            elif entity2_is_ncRNA:
                tax_ncRNA_protein_interactions.add(entity_key)
            elif entity2_is_miRNA:
                tax_miRNA_mRNA_interactions.add(entity_key)
        elif entity1_is_ncRNA:
            if entity2_is_protein:
                tax_ncRNA_protein_interactions.add(entity_key)
            elif entity2_is_ncRNA:
                tax_ncRNA_ncRNA_interactions.add(entity_key)
            elif entity2_is_miRNA:
                tax_ncRNA_ncRNA_interactions.add(entity_key)
        elif entity1_is_miRNA:
            if entity2_is_protein:
                tax_miRNA_mRNA_interactions.add(entity_key)
            elif entity2_is_ncRNA:
                tax_ncRNA_ncRNA_interactions.add(entity_key)
            elif entity2_is_miRNA:
                raise ValueError('Found miRNA-miRNA interaction: {} and {}.'.format(entity1, entity2))

    if len(assumed_ncRNA) > 0:
        logger.warning('Following {:d} entities were not classified as miRNA/ncRNA/protein. '
                        'Assuming they are ncRNA: {}'.format(len(assumed_ncRNA), ';'.join(sorted(list(assumed_ncRNA)))))
    if len(ppi_interactions) > 0:
        logger.warning('Found {:d} protein-protein interactions: {}'.format(len(ppi_interactions),
                                                                             '; '.join(sorted(list(ppi_interactions)))))

    tax_to_miRNA_mRNA_interactions[tax] = tax_miRNA_mRNA_interactions
    tax_to_ncRNA_ncRNA_interactions[tax] = tax_ncRNA_ncRNA_interactions
    tax_to_ncRNA_protein_interactions[tax] = tax_ncRNA_protein_interactions

tax_interaction_count.sort(key=lambda y: y[1], reverse=True)

name_to_tax = stringrnautils.species_name_to_taxonomy_id()
tax_to_name = dict((v, k) for k, v in name_to_tax.iteritems())

allowed_organisms = set()
sys.stderr.write('TaxonomyID\tSpeciesName\tUniqueInteractions\n')
for tax, interaction_count in tax_interaction_count:
    if tax in tax_to_name:
        species = tax_to_name[tax]
    else:
        species = 'NA'
    if interaction_count >= min_interaction_count:
        sys.stdout.write('%s ' % tax)
        allowed_organisms.add(tax)
    sys.stderr.write('%s\t%s\t%d\n' % (tax, species, interaction_count))

sys.stderr.write('\n\n')
sys.stderr.write('& miRNA-mRNA & ncRNA-protein & ncRNA-ncRNA & Sum\\\\\n')
sys.stderr.write('\\midrule\n')
# write Tex-style table that can be put into the supplement
for tax in [tup[0] for tup in tax_interaction_count]:
    if tax not in allowed_organisms:
        continue
    if tax in tax_to_name:
        species = tax_to_name[tax]
    else:
        species = 'NA'
    row = (species.capitalize(), 
           str(len(tax_to_miRNA_mRNA_interactions[tax])),
           str(len(tax_to_ncRNA_protein_interactions[tax])),
           str(len(tax_to_ncRNA_ncRNA_interactions[tax])),
           str(len(tax_to_miRNA_mRNA_interactions[tax]) +
               len(tax_to_ncRNA_protein_interactions[tax]) +
               len(tax_to_ncRNA_ncRNA_interactions[tax]))
           )
    sys.stderr.write(' & '.join(row))
    sys.stderr.write('\\\\\n')
sys.stderr.write('\\midrule\n')
row = ['Sum']
total_total = 0
for m in (tax_to_miRNA_mRNA_interactions, tax_to_ncRNA_protein_interactions, tax_to_ncRNA_ncRNA_interactions):
    int_count = 0
    for k, v in m.iteritems():
        if k in allowed_organisms:
            int_count += len(v)

    row.append(str(int_count))
    total_total += int_count
row.append(str(total_total))
sys.stderr.write(' & '.join(row))
sys.stderr.write('\\\\\n')

sys.stderr.write('\n\n')
sys.stderr.write('Organisms & miRNAs & ncRNAs (miRNAs excluded) & mRNAs/proteins & Sum\\\\\n')
sys.stderr.write('\\midrule\n')
# write Tex-style table that can be put into the supplement
for tax in [tup[0] for tup in tax_interaction_count]:
    if tax not in allowed_organisms:
        continue
    if tax in tax_to_name:
        species = tax_to_name[tax]
    else:
        species = 'NA'
    row = [species.capitalize()]
    total_count = 0
    for int_count in (len(tax_to_miRNAs[tax]), len(tax_to_ncRNAs[tax]), len(tax_to_proteins[tax])):
        row.append(str(int_count))
        total_count += int_count
    row.append(str(total_count))
    sys.stderr.write(' & '.join(row))
    sys.stderr.write('\\\\\n')

sys.stderr.write('\\midrule\n')
row = ['Sum']
total_total = 0
for m in (tax_to_miRNAs, tax_to_ncRNAs, tax_to_proteins):
    
    int_count = 0
    for k, v in m.iteritems():
        if k in allowed_organisms:
            int_count += len(v)
            
    row.append(str(int_count))
    total_total += int_count
row.append(str(total_total))
sys.stderr.write(' & '.join(row))
sys.stderr.write('\\\\\n')

evidence_channels = list(tax_to_source_to_interactions['9606'].keys())
evidence_channels.sort()
evidence_channel_table_names = [evidence_channel_to_table_name[x] for x in evidence_channels]

source_to_interaction_count = collections.defaultdict(int)
sys.stderr.write('\n\n')
sys.stderr.write('Organisms & ' + ' & '. join(evidence_channel_table_names) + ' & Sum\\\\\n')
sys.stderr.write('\\midrule\n')
# write Tex-style table that can be put into the supplement
for tax in [tup[0] for tup in tax_interaction_count]:
    if tax not in allowed_organisms:
        continue
    if tax in tax_to_name:
        species = tax_to_name[tax]
    else:
        species = 'NA'
    row = [species.capitalize()]

    row_total_set = set()
    for source in evidence_channels:
        if source in tax_to_source_to_interactions[tax]:
            int_count = len(tax_to_source_to_interactions[tax][source])
            row_total_set.update(tax_to_source_to_interactions[tax][source])
        else:
            int_count = 0
        row.append(str(int_count))
        source_to_interaction_count[source] += int_count
    row.append(str(len(row_total_set)))
    sys.stderr.write(' & '.join(row))
    sys.stderr.write('\\\\\n')
sys.stderr.write('\\midrule\n')
row = ['Sum']
total_total = 0
for source in evidence_channels:
    int_count = source_to_interaction_count[source]
    row.append(str(int_count))
    total_total += int_count
# row.append(str(total_total))
row.append(' ')
sys.stderr.write(' & '.join(row))
sys.stderr.write('\\\\\n')
logger.info('Done.' + os.linesep + os.linesep)
