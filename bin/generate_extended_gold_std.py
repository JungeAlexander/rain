#!/usr/bin/env python
import argparse
import os
import operator
from stringrnautils import Interaction, get_string_to_alias_mapper, get_unique_mir_mapper, EntityType, InteractionType, \
    get_non_coding_rna_alias_mapper
from collections import defaultdict, Counter
import logging

logger = logging.getLogger(os.path.basename(__file__))
logging.basicConfig(level=logging.INFO)

__author__ = 'Alexander Junge (alexander.junge@gmail.com)'


def is_existing_file(arg):
    if not os.path.isfile(arg):
        parser.error("Given file %s does not exist. Please provide an existing file name." % arg)
    return arg

parser = argparse.ArgumentParser(description="""
                                  Creates an extended gold standard master file.

                                  The extended gold standard contains interactions from the Croft et al. set
                                  supplemented with high confidence interactions from miRTarBase+NPInter combination.
                                  A high confidence interaction is an interaction that is an interaction supported
                                  by at least one low-throughput study and two additional studies.
                                  A low-throughput study is a study contributing at maximum 5 interactions to the
                                  miRTarBase+NPInter interaction set.
                                  """)
parser.add_argument('croft_master_file', type=is_existing_file,
                    help='The input master file derived from the Croft et al. set.')
parser.add_argument('miRTarBase_NPInter_master_file', type=is_existing_file,
                    help='The input master file derived from combining miRTarBase and NPInter.')
parser.add_argument('--human_only', action='store_true',
                    help='if set, only human interactions will be added to the extended gold standard')
parser.add_argument('--low-throughput-threshold', type=int,
                    help='a lowthroughput study is one with less interactions than this', default=5)
parser.add_argument('--n-low-throughput', type=int,
                    help='minimum number of lowtroughput experiments needed', default=1)
parser.add_argument('--n-total', type=int,
                    help='minimum total number of studies needed needed', default=3)
parser.add_argument('--low-throughput-pids', type=str,
                    help='a comma separated list of "high confidence" pids')


args = parser.parse_args()
n_requred = int(args.n_total)
low_throughput_threshold = int(args.low_throughput_threshold)
n_low_throughput_requred = int(args.n_low_throughput)

pid_file = open(args.low_throughput_pids, 'w')

logger.info("threashold = {0} n_low = {1} n_total = {2}".format(low_throughput_threshold, n_low_throughput_requred, n_requred))

logger.info('Generating extended gold standard (Croft + miRTarBase/NPInter high confidence).')

# map master files to interactions
in_master_files = [args.croft_master_file, args.miRTarBase_NPInter_master_file]
in_master_file_to_interactions = defaultdict(list)
pmid_to_miRTarBase_NPInter_interactions = defaultdict(int)
for master_file in in_master_files:
    logger.info('Reading master file %s.' % master_file)
    with open(master_file, 'r') as mf:
        for line in mf:
            org, ent1, ent2, directed, channel, score, sources, url, comment = line.rstrip('\n\r').split('\t')
            interact = Interaction(org, ent1, ent2, directed, channel, score, sources, url, comment)
            in_master_file_to_interactions[master_file].append(interact)
            if master_file == args.miRTarBase_NPInter_master_file:
                if url.startswith('http://www.ncbi.nlm.nih.gov/pubmed/'):
                    pubmed_ids = url[len('http://www.ncbi.nlm.nih.gov/pubmed/'):].split(',')
                    for pmid in pubmed_ids:
                        pmid_to_miRTarBase_NPInter_interactions[pmid] += 1
    logger.info('Found %d interactions in master file %s.' % (len(in_master_file_to_interactions[master_file]),
                                                               master_file))

# classify PubMed IDs as low throughput study
low_throughput_study_pmids = set()
for pmid, occ in pmid_to_miRTarBase_NPInter_interactions.items():
    if occ <= low_throughput_threshold:
        low_throughput_study_pmids.add(pmid)
logger.info('%d studies were classified as low throughput studies in miRTarBase/NPInter.' %
             len(low_throughput_study_pmids))
pid_file.write('{0}\n'.format(','.join(low_throughput_study_pmids)))

# collect all output interactions in a list; first add Croft et al.
output_interactions = in_master_file_to_interactions[args.croft_master_file]

# find high confidence interactions in miRTarBase/NPInter set of interactions
high_confidence_interaction_count = 0
for intact in in_master_file_to_interactions[args.miRTarBase_NPInter_master_file]:
    if intact._url.startswith('http://www.ncbi.nlm.nih.gov/pubmed/'):
        pubmed_ids = intact._url[len('http://www.ncbi.nlm.nih.gov/pubmed/'):].split(',')
        n_low_throughput_support = sum([x in low_throughput_study_pmids for x in pubmed_ids])
        if n_low_throughput_support >= n_low_throughput_requred and len(pubmed_ids) >= n_requred:
            output_interactions.append(intact)
            high_confidence_interaction_count += 1
logger.info('%d high confidence interactions were found in miRTarBase/NPInter.' % high_confidence_interaction_count)

# make sure that all interactions are unique using the fact that Interaction class implements __eq__ and __hash__ method
seen = set()
seen_add = seen.add
gold_standard_interactions = [x for x in output_interactions if not (x in seen or seen_add(x))]
logger.info('Extended gold standard contains %d interactions after redundancy filtering (before: %d)' %
             (len(gold_standard_interactions), len(output_interactions)))

if args.human_only:
    # reduce to human interactions
    logger.info('Filtering gold standatd to human interactions only.')
    human_gold_standard_interactions = [intact for intact in gold_standard_interactions
                                        if intact._org.startswith('9606')]
    gold_standard_interactions = human_gold_standard_interactions
    logger.info('Extended gold standard contains %d human interactions.' % len(human_gold_standard_interactions))

# print taxonomy distribution
taxonomy_id_to_occurrences = Counter([intact._org for intact in gold_standard_interactions])
sorted_taxonomy_ids_occurrences = sorted(taxonomy_id_to_occurrences.items(), key=operator.itemgetter(1), reverse=True)
logger.info('Taxonomy identifier\tOccurrences')
logger.info('--------------------------------')
for tax_id, occ in sorted_taxonomy_ids_occurrences:
    logger.info('%s\t%d' % (tax_id, occ))

# alias mappers used to determine
protein_alias_mapper = get_string_to_alias_mapper(taxonomy_id_to_occurrences.keys(), '', '', 10, 'all', True)
mir_alias_mapper = get_unique_mir_mapper()
ncrna_alias_mapper = get_non_coding_rna_alias_mapper()


def get_entity_type(entity_name, taxonomy_id):
    tax_protein_alias_mapper = protein_alias_mapper[taxonomy_id]
    entity_is_ncrna = entity_name in ncrna_alias_mapper[taxonomy_id]
    if entity_is_ncrna:
        return EntityType.ncRNA
    entity_is_protein = entity_name in tax_protein_alias_mapper
    if entity_is_protein:
        return EntityType.Protein
    entity_is_mirna = entity_name in mir_alias_mapper
    if entity_is_mirna:
        return EntityType.miRNA


# determine entity types for each interaction and print distribution
interaction_types = []
for intact in gold_standard_interactions:
    ent_1_type = get_entity_type(intact._ent1, intact._org)
    ent_2_type = get_entity_type(intact._ent2, intact._org)
    interaction_type = InteractionType.entities_to_interaction_type(ent_1_type, ent_2_type)
    interaction_types.append(interaction_type)
interaction_types_to_occurences = Counter(interaction_types)
sorted_interaction_type_occurrences = sorted(interaction_types_to_occurences.items(), key=operator.itemgetter(1),
                                             reverse=True)
logger.info('Interaction Type\tOccurrences')
logger.info('--------------------------------')
for int_type, occ in sorted_interaction_type_occurrences:
    logger.info('%s\t%d' % (InteractionType.interaction_type_to_string(int_type), occ))

# write extended gold standard to standard output
for intact in gold_standard_interactions:
    print(intact)
logger.info('Done.' + os.linesep + os.linesep)
