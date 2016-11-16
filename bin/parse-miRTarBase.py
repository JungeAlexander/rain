#!/usr/bin/env python
import argparse
import re
import stringrnautils
import os
import xlrd  # third party, can be installed by invoking: pip install xlrd
import logging

logger = logging.getLogger(os.path.basename(__file__))
logging.basicConfig(level=logging.INFO)

__author__ = 'Alexander Junge, Xiaoyong Pan'


def is_not_valid_path(arg):
    if os.path.exists(arg):
        parser.error("Given path %s exists. Please provide a non-existing output path or delete the existing file %s "
                     "manually before running this script." % (arg, arg))
    return arg


def is_valid_path(arg):
    if not os.path.exists(arg):
        parser.error("Given master file path %s does not exist." % arg)
    return arg


parser = argparse.ArgumentParser(description="""
                                 Parses miRNA-mRNA interactions from miRTarBase.
    All mRNAs are mapped to ENSPs using the STRING alias file whenever a mapping can be found, if no mapping is
    found, the corresponding interaction is ignored.
    miRNA identifiers are mapped to those miRBase identifiers used in RAIN.
    All identifier mapping is performed using the stringrnautils module.
    Interactions of miRNAs and mRNAs from different species are filtered out.
    All types of experiments supporting a certain interaction are first mapped to a 'cleaned' experiment name to
    remove typos or alternative names for the same experimental approach.
    """)
parser.add_argument('output_file_path', type=is_not_valid_path,
                    help='The output miRTarBase master file is written to this path.')
parser.add_argument('gold_standard_file', type=is_valid_path,
                    help='The gold standard file to benchmark against.')
parser.add_argument('-e', '--exclude_starbase',
                    help='Excludes interactions with PubMed IDs coming from StarBase from the output master file.',
                    action='store_true')
parser.add_argument('-i', '--include_starbase',
                    help='Includes interactions with PubMed IDs coming from StarBase in the output master file.',
                    action='store_true')
parser.add_argument('-d', '--data_path',
                    help='Path to data directory',
                    default="data")

args = parser.parse_args()

output_file_path = args.output_file_path
include_starbase = args.include_starbase
exclude_starbase = args.exclude_starbase
gold_standard_file_path = args.gold_standard_file
DATA_PATH = args.data_path

if (not include_starbase and not exclude_starbase) or (include_starbase and exclude_starbase):
    raise IOError("Please select exactly one of the following parameters: --exclude_starbase or --include_starbase.")
logger.info("Parsing miRTarBase.")

# ## CONFIGURATION

COMMENT = ''

# The scheme used to score interactions. Can be any of the following: 'PubMedIDs', 'assays' or 'evidence'. If None is
# given,no scoring is performed and interactions get score 0.
scoring_scheme = None

# Download miRTarBase
miRTarBase_version = "6.1"
if miRTarBase_version == "4.5":
    miRTarBase_file = os.path.join(DATA_PATH, "miRTarBase_" + miRTarBase_version + "_MTI.xls")
elif miRTarBase_version == "6.1":
    miRTarBase_file = os.path.join(DATA_PATH, "miRTarBase_" + miRTarBase_version + "_MTI.xlsx")
else:
    raise ValueError("miRTarBase version %s is not supported." % miRTarBase_version)

# INPUT FILES
# miRTarBaseFileName = 'data/miRTarBase_MTI.csv'  # miRTarBase xls file converted to csv file
# assayMappingFile = 'data/miRTarBase_assay_mapping.tsv'  # Maps assay names to 'cleaned' assay names
scoreMappingFile = os.path.join(DATA_PATH, "miRTarBase_assay_quality.tsv")  # Maps 'cleaned' assay names to conf. scores

# Maps species names in miRTarBase to corresponding taxonomy identifiers
uniqueSpeciesMap = {'Bombyx mori': '7091',
                    'Arabidopsis thaliana': '3702',
                    'Bos taurus': '9913',
                    'Caenorhabditis elegans': '6239',
                    'Xenopus laevis': '8355',
                    'Danio rerio': '7955',
                    'Drosophila melanogaster': '7227',
                    'Human cytomegalovirus': '10359',
                    'Homo sapiens': '9606',
                    'Gallus gallus': '9031',
                    'Kaposi sarcoma-associated herpesvirus': '37296',
                    'Mus musculus': '10090',
                    'Oryzias latipes': '8090',
                    'Ovis aries': '9940',
                    'Rattus norvegicus': '10116',
                    'Xenopus tropicalis': '8364',
                    'Epstein Barr virus': '10376',
                    'Mareks disease virus': '10390',
                    'Canis familiaris': '9615',
                    'Macaca mulatta': '9544',
                    'Solanum lycopersicum': '4081',
                    'Sus scrofa': '9832',
                    'Taeniopygia guttata': '59729'
                    }


def get_assay_mapping(assay_mapping_path):
    assay_dict = {}
    with open(assay_mapping_path, 'r') as fin:
        for curr_line in fin:
            split_cols = curr_line.rstrip().split('\t')
            assay_dict[split_cols[0]] = split_cols[1]
    return assay_dict


miRNA2Clean = stringrnautils.get_unique_mir_mapper()
miRNA2taxonomyID = stringrnautils.get_mir_id_to_tax_id_mapper()
targetName2targetID = stringrnautils.get_alias_to_string_mapper(organisms=uniqueSpeciesMap.values(),
                                                                filter_string_alias='', filter_string_id='')
restricted_pmids = stringrnautils.starbase_exp_pmids()

assayMappingFile = os.path.join(DATA_PATH, 'miRTarBase_assay_mapping.tsv')  # Maps assay names to 'cleaned' assay names
assay2Clean = get_assay_mapping(assayMappingFile)

not_mapped = 0
totalCount = 0

# Maps a certain interaction to a set of experiments supporting that interaction, the PubMedIDs and the evidence levels
# according to miRTarBase
interaction_to_experiments = {}
interaction_to_pubmed_ids = {}
interaction_to_evidence_types = {}

miRTarBase_workbook = xlrd.open_workbook(miRTarBase_file)
miRTarBase_sheet = miRTarBase_workbook.sheet_by_index(0)
num_rows = miRTarBase_sheet.nrows - 1
num_cells = miRTarBase_sheet.ncols - 1
curr_row = -1

cross_species_interactions = []
unknown_species = set()
not_mapped_mirnas = set()
wrong_tax_mirnas = set()
not_clean_assays = set()
interactions_in_starbase = 0
unmapped_mrnas = set()
not_string_mrnas = set()

while curr_row < num_rows:
    curr_row += 1
    # Skip header
    if curr_row == 0:
        continue
    row = miRTarBase_sheet.row(curr_row)
    #  File contains following columns:
    #  miRTarBase ID, miRNA, Species (miRNA), Target Gene, Target Gene (Entrez Gene ID), Species (Target Gene),
    #  Experiments, Support Type, References (PMID)
    cols = []
    curr_cell = -1
    while curr_cell < num_cells:
        curr_cell += 1
        try:
            str_cell = str(miRTarBase_sheet.cell_value(curr_row, curr_cell))
        except UnicodeEncodeError:
            str_cell = miRTarBase_sheet.cell_value(curr_row, curr_cell).encode('ascii', 'replace')
            logger.warning('Encoding entry in row %d cell %d failed. Using following representation: %s' %
                           (curr_row, curr_cell, str_cell))
        cols.append(str_cell)
    # PubMedID in last column is interpreted as float, so remove '.0' part
    cols[-1] = cols[-1][:-2]

    totalCount += 1
    if cols[2] != cols[5]:
        not_mapped += 1
        cross_species_interactions.append(cols[2] + ' and ' + cols[5])
        continue

    if cols[2] in uniqueSpeciesMap:
        organism = uniqueSpeciesMap[cols[2]]
    else:
        not_mapped += 1
        unknown_species.add(cols[2])
        continue

    # Map miRNA identifier to new miRBase identifiers
    if cols[1] in miRNA2Clean:
        miRNAId = miRNA2Clean[cols[1]]
    else:
        not_mapped += 1
        not_mapped_mirnas.add(cols[1])
        continue

    if miRNA2taxonomyID[miRNAId] != organism:
        wrong_tax_mirnas.add(miRNAId + ' is classified as ' + organism + ' but is ' + miRNA2taxonomyID[miRNAId])
        continue

    # Split into separate assays, map them to clean names
    assaysSplit = re.split("[,;] *|//", cols[6])  # 's/[,;] */\n/g; s/\/\//\n/g;'
    sep = ";"
    assays_clean = []
    idx = 0

    for ass in assaysSplit:
        if ass not in assay2Clean:
            not_clean_assays.add(ass)
        else:
            clean = assay2Clean[ass]
            assays_clean.append(clean)

    pubmed_id = cols[-1]
    if exclude_starbase:
        if int(pubmed_id) in restricted_pmids:
            not_mapped += 1
            interactions_in_starbase += 1
            continue
    elif include_starbase:
        pass
    else:
        raise IOError("--include_starbase or --exclude_starbase must be specified.")

    evidence_type = cols[-2].lower()

    mrna_name = cols[3]
    if organism in targetName2targetID:
        organism_map = targetName2targetID[organism]
        if mrna_name in organism_map:
            target_identifier = organism_map[mrna_name]
            triple = organism + '|' + miRNAId + '|' + target_identifier
            if triple in interaction_to_experiments:
                interaction_to_experiments[triple].update(assays_clean)
                interaction_to_pubmed_ids[triple].update([pubmed_id])
                interaction_to_evidence_types[triple].update([evidence_type])
            else:
                interaction_to_experiments[triple] = set(assays_clean)
                interaction_to_pubmed_ids[triple] = {pubmed_id}
                interaction_to_evidence_types[triple] = {evidence_type}
        else:
            not_mapped += 1
            unmapped_mrnas.add(mrna_name + '(' + organism + ')')
            continue
    else:
        not_mapped += 1
        not_string_mrnas.add(mrna_name + '(' + organism + ')')
        continue

# Determine scores using stringrnautils' benchmarking function
organisms = []
miRNAs = []
targets = []
experiments = []
pmid_bins = []
evidence_levels = []
for triple in interaction_to_experiments.keys():
    organism, miRNAId, targetIdentifier = triple.split('|')
    organisms.append(organism)
    miRNAs.append(miRNAId)
    targets.append(targetIdentifier)
    experiments.append(interaction_to_experiments[triple])
    pmid_count = len(interaction_to_pubmed_ids[triple])
    if pmid_count > 1:
        pmid_bin = "multiple PubMed IDs"
    else:
        pmid_bin = "single PubMed ID"
    # pmid_bin = str(pmid_count)
    pmid_bins.append(pmid_bin)
    evidence_types = interaction_to_evidence_types[triple]
    weak_evidence_types = [x for x in evidence_types if 'weak' in x]
    if len(weak_evidence_types) == len(evidence_types):
        evidence_levels.append('weak')
    else:
        evidence_levels.append('strong')

if not scoring_scheme:
    confidence_scores = [0] * len(organisms)
elif scoring_scheme == "PubMedIDs":
    confidence_scores = stringrnautils.benchmark(organisms, miRNAs, targets, pmid_bins, gold_standard_file_path,
                                                 discrete=True, fit_name='mirtarbase_pmids')
elif scoring_scheme == "assays":
    confidence_scores = stringrnautils.benchmark(organisms, miRNAs, targets, experiments, gold_standard_file_path,
                                                 discrete=True, fit_name='mirtarbase_assays')
elif scoring_scheme == "evidence":
    confidence_scores = stringrnautils.benchmark(organisms, miRNAs, targets, evidence_levels, gold_standard_file_path,
                                                 discrete=True, fit_name='mirtarbase_evidence')
else:
    raise IOError("Scoring scheme " + str(scoring_scheme) + " was not recognized.")

f = open(output_file_path, 'w')
for org, miR, tar, conf in zip(organisms, miRNAs, targets, confidence_scores):
    trip = org + '|' + miR + '|' + tar
    pubmed_url = 'http://www.ncbi.nlm.nih.gov/pubmed/' + ','.join(interaction_to_pubmed_ids[trip])
    # Write interaction
    f.write(org + '\t' + miR + '\t' + tar + '\t' + '1' + '\t' + 'Experiments' +
            '\t' + str(conf) +
            '\t' + 'miRTarBase' + '\t' + pubmed_url + '\t' + COMMENT + '\n')
f.close()

if len(cross_species_interactions) > 0:
    logger.info('Found {:d} cross-species interactions.'.format(len(cross_species_interactions)))
    logger.debug('Found {:d} cross-species interactions: {}'.format(len(cross_species_interactions),
                                                                    '; '.join(cross_species_interactions)))

if len(unknown_species) > 0:
    logger.info('Encountered {:d} unknown species identifiers.'.format(len(unknown_species)))
    logger.debug('Encountered {:d} unknown species identifiers: {}'.format(len(unknown_species),
                                                                           ', '.join(sorted(list(unknown_species)))))

if len(not_mapped_mirnas) > 0:
    logger.info('{:d} miRNA identifiers could not be mapped to miRBase aliases.'.format(len(not_mapped_mirnas)))
    logger.debug('{:d} miRNA identifiers could not be mapped to miRBase aliases: {}'.format(len(not_mapped_mirnas),
                                                                                            ', '.join(
                                                                                                sorted(list(
                                                                                                    not_mapped_mirnas)))))
if len(wrong_tax_mirnas) > 0:
    logger.info('Found {:d} miRNAs with wrong taxonomy ID.'.format(len(wrong_tax_mirnas)))
    logger.debug('Found {:d} miRNAs with wrong taxonomy ID: {}'.format(len(wrong_tax_mirnas),
                                                                     '; '.join(sorted(list(wrong_tax_mirnas)))))

if len(not_clean_assays) > 0:
    logger.info('No clean name found for {:d} assays.'.format(len(not_clean_assays)))
    logger.debug('No clean name found for {:d} assays: {}'.format(len(not_clean_assays),
                                                                    ', '.join(sorted(list(not_clean_assays)))))

if len(unmapped_mrnas) > 0:
    logger.info('Was not able to map {:d} mRNAs to STRING IDs.'.format(len(unmapped_mrnas)))
    logger.debug('Was not able to map {:d} mRNAs to STRING IDs: {}'.format(len(unmapped_mrnas),
                                                                            '; '.join(sorted(list(unmapped_mrnas)))))

if len(not_string_mrnas) > 0:
    logger.info('Was not able to map {:d} mRNAs to STRING IDs as STRING does not over organism.'.format(
        len(unmapped_mrnas)))
    logger.debug('Was not able to map {:d} mRNAs to STRING IDs as STRING does not over organism: {}'.format(
        len(unmapped_mrnas), '; '.join(sorted(list(unmapped_mrnas)))))

logger.info('Not added to RAIN: ' + str(not_mapped) + ' out of ' + str(totalCount) + ' interactions in miRTarBase')
if interactions_in_starbase > 0:
    logger.info('Of those, {:d} interactions were removed since their PubMed IDs were present in StarBase 2.0.'.format(
        interactions_in_starbase))

logger.info('Done.' + os.linesep + os.linesep)
