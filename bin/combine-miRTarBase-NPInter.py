#!/usr/bin/env python
import argparse
import stringrnautils
import collections
import os
import logging

logger = logging.getLogger(os.path.basename(__file__))
logging.basicConfig(level=logging.INFO)

__author__ = 'Alexander Junge (alexander.junge@gmail.com)'


def is_valid_path(arg):
    if not os.path.exists(arg):
        parser.error("Given master file path %s does not exist." % arg)
    return arg


def is_non_existing_path(arg):
    if os.path.exists(arg):
        parser.error("Given path %s exists. Please provide a non-existing output path or delete the existing file %s "
                     "manually before running this script." % (arg, arg))
    return arg


def is_existing_file(arg):
    if not os.path.isfile(arg):
        parser.error("Given file %s does not exist. Please provide an existing file name." % arg)
    return arg


parser = argparse.ArgumentParser(description=
                                 """Merges the miRTarBase and NPInter master files in one common set of interactions.
                                 Discrete scoring is performed where the bin assigned to each interaction is the number
                                 of PubMed IDs supporting the interaction.""")
parser.add_argument('miRTarBase_master_file', type=is_existing_file,
                    help='The input miRTarBase master file.')
parser.add_argument('NPInter_master_file', type=is_existing_file,
                    help='The input NPInter master file.')
parser.add_argument('output_master_file', type=is_non_existing_path,
                    help='The generated miRTarBase/NPInter master file is written to this path.')
parser.add_argument('gold_standard_file', type=is_valid_path,
                    help='The gold standard file to benchmark against.')
args = parser.parse_args()

logger.info("Combining miRTarBase and NPInter.")

##### CONFIGURATION
COMMENT = ''

miRTarBase_master_file = args.miRTarBase_master_file
NPInter_master_file = args.NPInter_master_file
out_master_file = args.output_master_file
gold_standard_file_path = args.gold_standard_file

min_pmid_count_in_last_bin = 3

master_files = [miRTarBase_master_file, NPInter_master_file]
database_names = ['miRTarBase', 'NPInter']

# Map PubMed IDs to number of interactions extracted from the given paper
pmid_occurrences = collections.defaultdict(int)
low_throughput_threshold = 5
medium_throughput_threshold = 30

interaction_to_sources = {}
interaction_to_pmids = {}
for index in range(len(master_files)):
    file_name = master_files[index]
    database_name = database_names[index]
    with open(file_name, 'r') as in_file:
        for line in in_file:
            cols = line.rstrip("\n").split('\t')
            if len(cols) != 9:
                raise IOError('Found wrong number of columns in master file ' + file_name)
            curr_key = (cols[0], cols[1], cols[2])
            curr_source = database_name
            curr_PubMedIDs = cols[-2].split('/')[-1].split(',')
            if curr_key in interaction_to_sources:
                interaction_to_sources[curr_key].update([curr_source])
                interaction_to_pmids[curr_key].update(curr_PubMedIDs)
            else:
                interaction_to_sources[curr_key] = {curr_source}
                interaction_to_pmids[curr_key] = set(curr_PubMedIDs)

orgs = []
id1s = []
id2s = []
pmid_bins = []
for triple in interaction_to_sources.keys():
    org, id1, id2 = triple
    orgs.append(org)
    id1s.append(id1)
    id2s.append(id2)
    pmid_count = len(interaction_to_pmids[triple])

    for pmid in interaction_to_pmids[triple]:
        pmid_occurrences[pmid] += 1

    if pmid_count < min_pmid_count_in_last_bin:
        pmid_bin = str(pmid_count)
    else:
        pmid_bin = str(min_pmid_count_in_last_bin) + '+'
    pmid_bins.append(pmid_bin)

scores = stringrnautils.benchmark(orgs, id1s, id2s, pmid_bins, gold_standard_file_path, discrete=True)

bin2score = {}
f = open(out_master_file, 'w')
for org, miR, tar, conf, curr_bin in zip(orgs, id1s, id2s, scores, pmid_bins):
    if not curr_bin in bin2score:
        bin2score[curr_bin] = conf
    trip = (org, miR, tar)
    sources = ','.join(interaction_to_sources[trip])
    pubmed_url = 'http://www.ncbi.nlm.nih.gov/pubmed/' + ','.join(interaction_to_pmids[trip])
    # Write interaction
    f.write(org + '\t' + miR + '\t' + tar + '\t' + '1' + '\t' + 'experimental' +
                  '\t' + str(conf) +
                  '\t' + sources + '\t' + pubmed_url + '\t' + COMMENT + '\n')
f.close()

# Remove input master files
for master_file in master_files:
    os.unlink(master_file)

# Classify PubMed IDs as low, medium and high throughput papers
low_alias = 3
medium_alias = 2
high_alias = 1
pmid_to_type = {}
for pmid, occ in pmid_occurrences.items():
    if occ <= low_throughput_threshold:
        pmid_to_type[pmid] = low_alias
    elif occ <= medium_throughput_threshold:
        pmid_to_type[pmid] = medium_alias
    else:
        pmid_to_type[pmid] = high_alias

bin2type2count = {}
for org, miR, tar, curr_bin in zip(orgs, id1s, id2s, pmid_bins):
    trip = (org, miR, tar)
    pmids = interaction_to_pmids[trip]
    curr_type = 0
    for pmid in pmids:
        curr_pmid_type = pmid_to_type[pmid]
        curr_type = max(curr_pmid_type, curr_type)

    if curr_bin in bin2type2count:
        bin2type2count[curr_bin][curr_type] += 1
    else:
        bin2type2count[curr_bin] = collections.defaultdict(int)
        bin2type2count[curr_bin][curr_type] += 1

logger.info('\t&\t'.join(
    ['Bin', 'Supported by >=1 LT', 'No LT support, but >=1 MT', 'Only HT support', 'Total', 'Current Score']))
for curr_bin in sorted(bin2type2count.keys()):
    total = 0
    for val in bin2type2count[curr_bin].values():
        total += val
    lt_frac = bin2type2count[curr_bin][low_alias]/float(total)
    mt_frac = bin2type2count[curr_bin][medium_alias]/float(total)
    ht_frac = bin2type2count[curr_bin][high_alias]/float(total)

    logger.info(curr_bin + '\t&\t' +
        str(bin2type2count[curr_bin][low_alias]) + '(' + str(round(lt_frac, 3)) + ')' + '\t&\t' +
        str(bin2type2count[curr_bin][medium_alias]) + '(' + str(round(mt_frac, 3)) + ')' + '\t&\t' +
        str(bin2type2count[curr_bin][high_alias]) + '(' + str(round(ht_frac, 3)) + ')' + '\t&\t' +
        str(total) + '\t&\t' +
        str(bin2score[curr_bin]))
logger.info("\nLT = 'Low-throughput experiment'. Experiments are considered low-throughput if they are mentioned as "
                 "references for %d or less interactions." % low_throughput_threshold)
logger.info("MT = 'Medium-throughput experiment'. Experiments are considered medium-throughput if they are "
                 "mentioned as references for more than %d but less or equal than %d experiments."
                 % (low_throughput_threshold, medium_throughput_threshold))
logger.info("HT = 'High-throughput experiment'. Experiments are considered high-throughput if they are "
                 "mentioned as references for more than %d experiments."
                 % (medium_throughput_threshold))

logger.info('Done.' + os.linesep + os.linesep)

