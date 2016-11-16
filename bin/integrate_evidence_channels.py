#!/usr/bin/env python
import argparse
import collections
import gzip
import os
import re
import sys
import logging

import operator

logger = logging.getLogger(os.path.basename(__file__))
logging.basicConfig(level=logging.INFO)

__author__ = 'Alexander Junge (alexander.junge@gmail.com)'


def float_geq_0_leq_1(arg):
    try:
        float_arg = float(arg)
        if float_arg >= 0:
            if float_arg <= 1:
                return float_arg
            else:
                parser.error('Given float {:f} is not <= 1.'.format(float_arg))
        else:
            parser.error('Given float {:f} is not >= 0.'.format(float_arg))
    except ValueError:
        parser.error('Given value {} is not a float.'.format(arg))


parser = argparse.ArgumentParser(description='''
Computes integrated reliability scores for interactions and prints those interactions to stdout whose integrated scores
pass a given threshold. Input interactions are read from stdin.
''')
parser.add_argument('score_threshold', type=float_geq_0_leq_1, help='Integrated reliability score threshold to use.')
parser.add_argument('--STRING_prior', type=float_geq_0_leq_1, default=0.04113021,
                    help='Prior probability used by STRING during score integration. Defaults to STRING 10 prior.')
parser.add_argument('--integrated_score_file_path', default=os.path.join('data', 'integrated_scores.tsv.gz'),
                    help='scores obtained after integration will be written to this file')
args = parser.parse_args()

score_threshold = args.score_threshold
string_prior = args.STRING_prior
integrated_score_file_path = args.integrated_score_file_path

logger.info('Computing integrated scores. Filtering interactions with integrated score < {:f}.'.format(score_threshold))
interactions_to_evidence_channel_to_score = collections.defaultdict(dict)
interactions_to_master_file_rows = collections.defaultdict(list)
alternative_keys_to_final_keys = {}

scores_below_prior = 0
for line in sys.stdin:
    if re.match('^\s*#', line):
        continue
    row = line.rstrip('\r\n').split('\t')

    try:
        organism, node1, node2, directed, evidence, score, source, link, comment = row
    except ValueError as e:
        logger.critical('Error (see below) encountered when unpacking row in aggregate master file. '
                        'Corresponding row was:\n%s' % row)
        raise e

    if float(score) < string_prior:  # ignore scores smaller than the prior as these would worsen the combined score
        scores_below_prior += 1
        continue

    key_1 = (node1, node2, organism)
    if key_1 in alternative_keys_to_final_keys:
        key = alternative_keys_to_final_keys[key_1]
    else:
        key = key_1
        key_2 = (node2, node1, organism)
        alternative_keys_to_final_keys[key_1] = key
        alternative_keys_to_final_keys[key_2] = key

    if evidence in interactions_to_evidence_channel_to_score[key]:
        msg = 'Interaction between {} and {} in organism {} appeared more than once in evidence channel {}. ' \
              'Exiting!'.format(key[0], key[1], key[2], evidence)
        logger.critical(msg)
        sys.exit(msg)

    interactions_to_evidence_channel_to_score[key][evidence] = float(score)
    interactions_to_master_file_rows[key].append(line)
logger.info('Ignores {:d} scores that were smaller than the prior {:f}.'.format(scores_below_prior, string_prior))

interactions_to_evidence_channel_to_integrated_score = {}
for key, evidence_dict in interactions_to_evidence_channel_to_score.iteritems():
    # (1 - p)   | |  (1 - Pi)
    # ------- = | |  --------
    # (1 - p*)  | |  (1 - Pi*)
    #
    # p = 1 - (1 - p*)^(1-N) * Product((1 - pi) --> this only applies when Pi*==p*
    rev_scores = [1 - s for s in evidence_dict.values()]
    rev_priors = [1 - string_prior] * len(evidence_dict)
    integrated_score = 1 - ((1 - string_prior) * reduce(operator.mul, rev_scores)/reduce(operator.mul, rev_priors))
    final_score = max(integrated_score, 0.0)  # ensure that min score is 0.0
    if final_score > 1:
        msg = 'Interaction between {} and {} in organism {} has integrated score > 1. ' \
              'Exiting!'.format(key[0], key[1], key[2], str(final_score))
        logger.critical(msg)
        sys.exit(msg)
    interactions_to_evidence_channel_to_integrated_score[key] = final_score

interactions_passing_threshold = 0
with gzip.open(integrated_score_file_path, 'w') as integrated_score_file:
    for key, integrated_score in interactions_to_evidence_channel_to_integrated_score.iteritems():
        integrated_score_file.write('\t'.join((key[2], ) + key[:2]) +
                                    '\t' + ';'.join(sorted(interactions_to_evidence_channel_to_score[key].iterkeys())) +
                                    '\t' + str(integrated_score) + '\n')
        if integrated_score >= score_threshold:
            for row in interactions_to_master_file_rows[key]:
                sys.stdout.write(row)
            interactions_passing_threshold += 1
logger.info('Integrated scores computed and filtered. {:d} out of {:d} interactions '
            'remain.'.format(interactions_passing_threshold,
                             len(interactions_to_evidence_channel_to_integrated_score)))
