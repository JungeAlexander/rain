"""
NAME
    parse-RAID.py

DESCRIPTION
    Parses RNA-Protein and RNA-RNA interactions from the RAID database, maps identifiers to ENSPs wherever possible and
    converts the database to master file format.
    Raw interaction scores in the master file are the number of literature references supporting a certain interaction.

AUTHOR
    Alexander Junge <alexander.junge@gmail.com>
"""

import sys
import argparse
import gzip
import stringrnautils

# Handle CLI arguments
parser = argparse.ArgumentParser(description='Parses RNA-Protein and RNA-RNA interactions from the RAID database, '
                                             'maps identifiers to ENSPs wherever possible and '
                                             'converts the database to master file format. '
                                             'Raw interaction scores in the master file are the number of literature '
                                             'references supporting a certain interaction.')
parser.add_argument('-p', '--rna-proteins', dest='parse_rna_proteins', help='Parse RAID RNA-protein interactions.',
                    action='store_true')
parser.add_argument('-r', '--rna-rna', dest='parse_rna_rna', help='Parse RAID RNA-protein interactions.',
                    action='store_true')
args = parser.parse_args()

# ## CONFIGURATION
# The RAID RNA-protein interaction file preprocessed by script RAID2ENSEMBLE_conversion.py
rna_protein_file_path = 'data/RAID-RNA-Protein-downloaded-01-August-2014.ENSG.tsv'
# The RAID RNA-RNA interaction file preprocessed by script RAID2ENSEMBLE_conversion.py
rna_rna_file_path = 'data/RAID-RNA-RNA-downloaded-01-August-2014.ENSG.tsv'

rna_protein_output_file_path = 'RAID-RNA-Protein-master-file.tsv'
rna_rna_output_file_path = 'RAID-RNA-RNA-master-file.tsv'

if args.parse_rna_proteins:
    input_file_path = rna_protein_file_path
    output_file_path = rna_protein_output_file_path
elif args.parse_rna_rna:
    input_file_path = rna_rna_file_path
    output_file_path = rna_rna_output_file_path
else:
    raise ValueError('Either -p (--rna-proteins) or -r (--rna-rna) switch must be set.')

string_aliases_file_path = 'data/protein.aliases.v9.1.txt.gz'
mirbase_aliases_file_path = 'data/mirBase_aliases20.tsv'
ensg_to_mirbase_file_path = 'data/EnsembleGeneID2miRBaseID.gz'

human_taxonomy_id = '9606'


def create_string_alias_to_ensp_dict(alias_file_path):
    """

    :param alias_file_path: a string: location of the alias file
    :return: dictionary: STRING alias name -> ENSP
    """
    sys.stderr.write('Started mapping gene names to ENSPs.\n')
    with gzip.open(alias_file_path, 'r') as alias_file:
        return_dict = {}
        for l in alias_file:
            if l[0] == '#':
                continue
            values = l.split('\t')
            if values[0] != human_taxonomy_id:
                continue
            key = values[2]
            return_dict[key] = values[1]
            self_key = values[1]
            if not self_key in return_dict:
                return_dict[self_key] = values[1]
        sys.stderr.write('Finished mapping gene names to ENSPs.\n')
        return return_dict


def create_ensg_to_mirbase_dict(ensg_mirbase_file_path):
    """

    :param ensg_mirbase_file_path: a string: The location of the ENSG to miRBase ID mapping file obtained from BioMart
    :return: a dictionary: ENSG identifier -> miRBase Identifier
    """
    mirbase_dict = {}
    with gzip.open(ensg_mirbase_file_path, 'r') as mirbase_file_handle:
        # Skip header
        mirbase_file_handle.readline()

        for l in mirbase_file_handle:
            cols = l.rstrip().split('\t')
            # Ignore those ENSGs that do not correspond to miRNAs
            if len(cols) < 2:
                continue
            ensg = cols[0]
            mirbase_id = cols[1]

            # Map miRBase ID to the one used in STRING RNA
            if mirbase_id in mirbase_aliases:
                mirbase_id = mirbase_aliases[mirbase_id]
                mirbase_dict[ensg] = mirbase_id
            else:
                sys.stderr.write('Was not able to obtain miRBase identifier for: ' + mirbase_id + '\n')
    return mirbase_dict


def get_final_identifier(identifier, molecule_type):
    """

    :param identifier: a string: identifier of the molecule
    :param molecule_type: a string: type of the molecule
    :return: the final identifier of the molecule as used in STRING RNA
    """
    if molecule_type == 'Protein' or molecule_type == 'mRNA':
        if identifier in string_aliases_to_ensp:
            new_identifier = string_aliases_to_ensp[identifier]
        else:
            sys.stderr.write(molecule_type + ' ' + identifier + ' not in STRING alias file.' + '\n')
            new_identifier = identifier
    elif molecule_type == 'other':
        # Allow if an ENSG is linked to the molecule
        if identifier.startswith('ENSG'):
            # There is only once case, which is a protein misclassified as other
            if identifier in string_aliases_to_ensp:
                new_identifier = string_aliases_to_ensp[identifier]
            else:
                new_identifier = identifier
        else:
            new_identifier = None
    elif molecule_type == 'lncRNA':
        if identifier.startswith('ENSG'):
            sys.stderr.write('\tFound lncRNA with identifier: ' + identifier + '\n')
        new_identifier = identifier
    elif molecule_type == 'snoRNA' or molecule_type == 'snRNA':
        new_identifier = identifier
    elif molecule_type == 'miRNA':
        if identifier in ensg_to_mirbase:
            new_identifier = ensg_to_mirbase[identifier]
        else:
            new_identifier = None
    else:
        sys.stderr.write('Unknown molecule type ' + molecule_type + ' with identifier ' + identifier + '\n')
        new_identifier = None
    return new_identifier


string_aliases_to_ensp = create_string_alias_to_ensp_dict(string_aliases_file_path)
mirbase_aliases = stringrnautils.get_mir_mapper(mirbase_aliases_file_path)
ensg_to_mirbase = create_ensg_to_mirbase_dict(ensg_to_mirbase_file_path)

# Maps each interaction of to a set of PubMed IDs supporting this interaction
interaction_to_pmids = {}
missing_interactions = 0
# First pass over file just gathers all PubMed IDs for all interactions
with open(input_file_path, 'r') as f:
    for line in f:
        columns = line.split('\t')
        first_id = columns[3]
        second_id = columns[7]

        if first_id == '' or second_id == '':
            missing_interactions += 1
            continue

        pmid = columns[11]
        interaction_key = first_id + second_id
        if interaction_key in interaction_to_pmids:
            interaction_to_pmids[interaction_key].add(pmid)
        else:
            new_set = set()
            new_set.add(pmid)
            interaction_to_pmids[interaction_key] = new_set
sys.stderr.write('Interactions omitted due to missing identifiers: ' + str(missing_interactions) + '\n')

# Contains all interactions that have already been written to the output master file
written_interactions = set()
interaction_count = 0
not_mapped_interactions = 0
with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
    for line in input_file:
        columns = line.split('\t')
        raid_id = columns[0]
        first_id = columns[3]
        first_type = columns[4]
        second_id = columns[7]
        second_type = columns[8]
        interaction_key = first_id + second_id

        # Make sure that each interaction is written only once
        if interaction_key in written_interactions:
            continue

        written_interactions.add(interaction_key)
        raw_score = len(interaction_to_pmids[interaction_key])
        interaction_count += 1

        pubmed_ids = ';'.join(interaction_to_pmids[interaction_key])

        # Convert to identifiers used in STRING RNA. If a mapping fails, omit this column.
        final_first_id = get_final_identifier(first_id, first_type)
        final_second_id = get_final_identifier(second_id, second_type)



        if final_first_id is None:
            sys.stderr.write(first_type + ' identifier ' + first_id + ' could not be mapped.\n')
        if final_second_id is None:
            sys.stderr.write(second_type + ' identifier ' + second_id + ' could not be mapped.\n')
        if final_first_id is None or final_second_id is None:
            not_mapped_interactions += 1
            continue

        # All fields present in the master file
        organism = human_taxonomy_id
        id1 = final_first_id
        id2 = final_second_id
        directed = '0'
        evidence = 'Databases'
        score = str(raw_score)
        source = 'RAID'
        url = 'http://www.rna-society.org/raid/detail.php?raid=' + raid_id

        # Write to the master file
        output_line = '\t'.join((organism, id1, id2, directed, evidence, score, source, url))
        # For printing the PubMed IDs:
        # output_line = id1 + '\t' + id2 + '\t' + pubmed_ids
        output_file.write(output_line + '\n')

sys.stderr.write('Number of unique interactions: ' + str(interaction_count) + '\n')
sys.stderr.write('Number of omitted interactions due to not mappable identifiers: ' + str(not_mapped_interactions) +
                 '\n')
