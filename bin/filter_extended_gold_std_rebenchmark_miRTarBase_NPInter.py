#!/usr/bin/env python
import os
from collections import defaultdict
import logging
import argparse
from stringrnautils import Interaction, benchmark

logger = logging.getLogger(os.path.basename(__file__))
logging.basicConfig(level=logging.INFO)

__author__ = 'Alexander Junge (alexander.junge@gmail.com)'


def is_existing_file(arg):
    if not os.path.isfile(arg):
        parser.error("Given file %s does not exist. Please provide an existing file name." % arg)
    return arg

parser = argparse.ArgumentParser(description="""
                                 Removes interactions from miRTarBase/NPInter master file that were added to extended
                                 gold standard interactions. Cleaned miRTarBase/NPInter master file is written to
                                 standard output. Also a re-benchmarking based on the extended gold standard is
                                 performed.
                                 """)
parser.add_argument('extended_gold_std_master_file', type=is_existing_file,
                    help='The input master file representing the extended gold standard set.')
parser.add_argument('miRTarBase_NPInter_master_file', type=is_existing_file,
                    help='The input master file derived from combining miRTarBase and NPInter.')
parser.add_argument('--min_pmid_count_in_last_bin', type=int, default=3,
                    help='The minimum number of PubMed IDs associated with an interaction to assign it to the bin with '
                         'highest reliability.')
parser.add_argument('--low-throughput-pid-file', type=is_existing_file,
                    help='path of the file where all low-throughput pids are stored')
parser.add_argument('--no-filter',  action='store_true', default=False)
args = parser.parse_args()

# import pdb;pdb.set_trace()
low_throughput_pids = set(open(args.low_throughput_pid_file).read().rstrip().split(','))

logger.info('Removing interactions from miRTarBase/NPInter that were added to extended gold standard interactions ' +
            'also re-benchmarks miRTarBase/NPInter.')
in_master_files = [args.extended_gold_std_master_file, args.miRTarBase_NPInter_master_file]
in_master_file_to_interactions = defaultdict(list)
for master_file in in_master_files:
    logger.info('Reading master file %s.' % master_file)
    with open(master_file, 'r') as mf:
        for line in mf:
            org, ent1, ent2, directed, channel, score, sources, url, comment = line.rstrip('\n\r').split('\t')
            interact = Interaction(org, ent1, ent2, directed, channel, score, sources, url, comment)
            in_master_file_to_interactions[master_file].append(interact)
    logger.info('Found %d interactions.' % len(in_master_file_to_interactions[master_file]))

# reduce extended gold standard to those interactions coming from miRTarBase/NPInter
gold_mitarbase_npinter = []
for intact in in_master_file_to_interactions[args.extended_gold_std_master_file]:
    source_lower = intact._sources.lower()
    if 'npinter' in source_lower or 'mirtarbase' in source_lower:
        gold_mitarbase_npinter.append(intact)
logger.info('Found %d interactions in gold standard that stem from miRTarBase/NPInter.' % len(gold_mitarbase_npinter))


################################################################################
# altered because we not do the Low throughput stuff
# therefore we keep all
new_mirtarbase_npinter_interactions = in_master_file_to_interactions[args.miRTarBase_NPInter_master_file]
################################################################################
# iterate over miRTarBase/NPInter interactions and keep only those that were not added to gold standard
# new_mirtarbase_npinter_interactions = []
# for intact in in_master_file_to_interactions[args.miRTarBase_NPInter_master_file]:
#     if intact not in gold_mitarbase_npinter:  # Interaction classs supports membership test in list of Interactions?!
#         new_mirtarbase_npinter_interactions.append(intact)
################################################################################

# re-benchmark miRTarBase/NPInter using extended gold standard
orgs = []
id1s = []
id2s = []
pmid_bins = []
high_throughput_pmid_bins=[]
for intact in new_mirtarbase_npinter_interactions:
    org = intact._org
    id1 = intact._ent1
    id2 = intact._ent2
    url = intact._url
    orgs.append(org)
    id1s.append(id1)
    id2s.append(id2)

    if url.startswith('http://www.ncbi.nlm.nih.gov/pubmed/'):
        pubmed_ids = url[len('http://www.ncbi.nlm.nih.gov/pubmed/'):].split(',')
        pmid_count = len(pubmed_ids)
        high_throughput_pmid_counts = len(set(pubmed_ids) - low_throughput_pids)
    else:
        raise ValueError('Found interaction in miRTarBase/NPInter that does not have a PubMed URI:\n%s' % str(intact))

    if pmid_count < args.min_pmid_count_in_last_bin:
        pmid_bin = str(pmid_count)
    else:
        pmid_bin = str(args.min_pmid_count_in_last_bin) + '+'

    if high_throughput_pmid_counts < args.min_pmid_count_in_last_bin:
        high_throughput_pmid_bin = str(high_throughput_pmid_counts)
    else:
        high_throughput_pmid_bin = str(args.min_pmid_count_in_last_bin) + '+'

    pmid_bins.append(pmid_bin)
    # if high_throughput_pmid_bin != "0":
    high_throughput_pmid_bins.append(high_throughput_pmid_bin)

if args.no_filter:
    high_throughput_pmid_bins = None

new_scores = benchmark(orgs, id1s, id2s, pmid_bins, args.extended_gold_std_master_file, discrete=True,
                       fit_name='rebenchmark_miRTarBase_NPInter', filtered_scores=high_throughput_pmid_bins)

for intact, new_score in zip(new_mirtarbase_npinter_interactions, new_scores):
    intact._score = str(new_score)

# print re-benchmarked miRTarBase/NPInter interactions
for intact in new_mirtarbase_npinter_interactions:
    print intact
logger.info('Done.' + os.linesep + os.linesep)