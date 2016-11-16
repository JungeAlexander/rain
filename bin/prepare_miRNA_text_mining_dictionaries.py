#!/usr/bin/env python
import argparse
import stringrnautils
import os
import re
import gzip
import logging

logger = logging.getLogger(os.path.basename(__file__))
logging.basicConfig(level=logging.INFO)

__author__ = 'Alexander Junge, RTH, University of Copenhagen <alexander.junge@gmail.com>'

parser = argparse.ArgumentParser(description="""
Prepares miRNA alias dictionaries and miRNA orthology information used by
RAIN's text mining pipeline.""")
args = parser.parse_args()

# Constants
MIRBASE_VERSION = stringrnautils.MIRBASE_VERSION
DICTIONARY_FOLDER = stringrnautils.ID_PATH
DATA_FOLDER = stringrnautils.DATA_PATH

# miRBase organism codes extracted from http://mirbase.org/help/genome_summary.shtml
MIR_BASE_SPECIES_PREFIXES = set()
MIR_BASE_SPECIES_PREFIXES.update(['aae', 'aca', 'aga', 'aly', 'ame', 'api', 'ata', 'ath', 'atr', 'bdi', 'bfl', 'blv',
                                  'bma', 'bmo', 'bna', 'bra', 'bta', 'cbn', 'cbr', 'cel', 'cfa', 'chi', 'cin', 'cqu',
                                  'cre', 'crm', 'csa', 'cte', 'dan', 'ddi', 'der', 'dgr', 'dme', 'dmo', 'dpe', 'dps',
                                  'dpu', 'dre', 'dse', 'dsi', 'dvi', 'dwi', 'dya', 'ebv', 'eca', 'efu', 'egr', 'emu',
                                  'esi', 'fru', 'gga', 'ggo', 'gma', 'gra', 'hbr', 'hbv', 'hcm', 'hme', 'hsa', 'isc',
                                  'ksh', 'lgi', 'lja', 'mdm', 'mdo', 'mdv', 'mes', 'mgh', 'mml', 'mmu', 'mse', 'mtr',
                                  'nve', 'nvi', 'oan', 'oar', 'ocu', 'oha', 'ola', 'osa', 'pma', 'ppc', 'ppe', 'ppt',
                                  'ppy', 'prd', 'ptc', 'ptr', 'pxy', 'rco', 'rno', 'rrv', 'sbi', 'sha', 'sja', 'sko',
                                  'sly', 'sma', 'sme', 'spu', 'ssc', 'str', 'stu', 'tca', 'tch', 'tgu', 'tni', 'vvi',
                                  'xtr', 'zma'])

# INPUT FILES
# location of miRBase family classification
MIR_FAMILY_FILE = os.path.join(DATA_FOLDER, 'miRBase_%s_miFam.dat.gz' % MIRBASE_VERSION)

# OUTPUT FILES
# Output file storing miRNA aliases used for text mining. Tab-delimited with three columns:
# taxonomy ID<tab>STRING RNA miRNA identifier<tab>STRING RNA miRNA alias
OUT_MIR_ALIAS_FILE = os.path.join(DICTIONARY_FOLDER, 'miRNA_text_mining_aliases_miRBase_%i.tsv.gz' % MIRBASE_VERSION)
# Output file mapping miRNA precursors to mature miRNAs generated from the precursor
OUT_MIR_PRECURSOR_TO_MATURE_FILE = os.path.join(DICTIONARY_FOLDER,
                                                'miRNA_precursor_to_mature_miRBase_%i.tsv.gz' % MIRBASE_VERSION)


def __create_mir_alias_dictionary__():
    """ Creates a miRNA alias dictionary file (at OUT_MIR_ALIAS_FILE) containing the following mappings:

        a) each RAIN miR identifier mapped to itself
        b) each miR's aliases mapped to the miR identifier used in RAIN
        c) the precursor miR mapped to the identifier of the derived mature miR used in RAIN if and only if only a
           single mature miR for the respective precursor is contained in miRBase
        d) all of the mappings in a), b), c) with the species prefix removed, e.g., hsa-miR-125a-5p becomes miR-125a-5p
    """
    # The miR alias dictionary to be filled. Maps each two-element tuple (taxonomy ID, alias) to a two element tuple
    # (taxonomy ID, RAIN identifier)
    # Example:
    # mir_alias_mapper[('9606','hsa-miR-125a*')] = ('9606', 'miR-125a-3p')
    mir_alias_mapper = {}

    # Dictionary used to obtain the mappings a) and b).
    unique_mir_mapper = stringrnautils.get_unique_mir_mapper()

    # Instead of simply using the unique_mir_mapper dictionary for the following,
    # the underlying mapping file is read to get access to the taxonomy IDs of miRs.
    default_mir_mapping_file_path = stringrnautils.MIR_MAPPING_ALIASES_PATH % MIRBASE_VERSION
    with open(default_mir_mapping_file_path, 'r') as default_mir_mapping_file:
        for line in default_mir_mapping_file:
            organism, target_mir, alias_mir, priority = line.rstrip().split('\t')
            # Following check only adds mappings for aliases that uniquely map to a miR identifier.
            # Note that this takes care of the condition c) required for the output mapping of this method.
            if alias_mir in unique_mir_mapper:
                mir_alias_mapper[(organism, alias_mir)] = (organism, target_mir)

    # Add the following to the output mapping:
    # d) all of the mappings in a), b), c) with the species prefix removed, e.g., hsa-miR-125a-5p becomes miR-125a-5p
    for taxon_alias, taxon_id_tuple in mir_alias_mapper.items():
        taxonomy_id, alias = taxon_alias
        alias_split = alias.split('-')
        if len(alias_split) > 1:
            alias_prefix = alias_split[0]
            if alias_prefix in MIR_BASE_SPECIES_PREFIXES:
                mir_alias_mapper[(taxonomy_id, '-'.join(alias_split[1:]))] = taxon_id_tuple

    # write the final miR mapping before returning
    with gzip.open(OUT_MIR_ALIAS_FILE, 'wb') as mir_alias_file:
        for taxon_alias_tuple, taxon_id_tuple in mir_alias_mapper.items():
            tax_alias, alias = taxon_alias_tuple
            taxonomy_id, mir_id = taxon_id_tuple
            mir_alias_file.write('%s\t%s\t%s\n' % (taxonomy_id, mir_id, alias))


def get_text_mining_mir_dictionary():
    """
    Generates and return miRNA aliases used for text mining. See __create_mir_alias_dictionary__ for a detailed
    description of this mapping.
    :return: a dictionary mapping a tuple of taxonomy identifiers and miRNA alias to final miRNA identifiers, e.g.,
             ('9606','hsa-miR-125a*') is mapped to 'miR-125a-3p'
    """
    if logger.getEffectiveLevel() == logging.DEBUG or not os.path.exists(OUT_MIR_ALIAS_FILE):
        __create_mir_alias_dictionary__()

    mir_alias_to_identifier = {}
    with gzip.open(OUT_MIR_ALIAS_FILE, 'rb') as mir_alias_file:
        for line in mir_alias_file:
            tax_id, mir_id, mir_alias = line.rstrip('\r\n').split('\t')
            mir_alias_to_identifier[(tax_id, mir_alias)] = mir_id
    return mir_alias_to_identifier


def __create_mir_precursor_to_mature_mapping_file__():
    """
    Creates a tab-delimited file at OUT_MIR_PRECURSOR_TO_MATURE_FILE with three columns:
    taxonomy ID<tab>miRBase precursor ID<tab>miRBase hairpin ID 1;miRBase hairpin ID 2; ...
    """
    organisms_found = set()  # for sanity checking, collect the names of all organisms found in the database
    organism_mapper = stringrnautils.species_name_to_taxonomy_id()  # extract mappings only for organisms in STRING

    organisms_not_mapped = set()
    with gzip.open(stringrnautils.MIRBASE_MIRNA_DAT_PATH % MIRBASE_VERSION, 'rb') as mirna_dat_file, \
            gzip.open(OUT_MIR_PRECURSOR_TO_MATURE_FILE, 'wb') as mir_precursor_to_mature_file:
        mature_accessions = []
        mature_identifiers = []
        precursor_accession = None
        precursor_identifier = None
        for line in mirna_dat_file:
            tag = line[:2]
            if tag == 'ID':  # the human readable identifier of the precursor
                precursor_identifier = line.split()[1]
                logger.debug('Found miRNA precursor %s.' % precursor_identifier)
            elif tag == 'AC':  # the accession of the precursor
                precursor_accession = line.split()[1].rstrip(';')
                logger.debug('Found miRNA precursor accession %s.' % precursor_accession)
            elif tag == 'DE':
                organism = ' '.join(line.split('stem')[0].split()[1:-1]).lower()
            elif tag == 'FT':  # information about mature miRNAs made from the current precursor
                mature_mir_accession_match = re.search(r'/accession=\"(\w+)\"', line)
                if mature_mir_accession_match:
                    mature_mir_accession = mature_mir_accession_match.groups()[0]
                    mature_accessions.append(mature_mir_accession)
                mature_mir_identifier_match = re.search(r'/product=\"(.+)\"', line)
                if mature_mir_identifier_match:
                    mature_identifiers.append(mature_mir_identifier_match.groups()[0])
            elif tag == '//':
                if organism in organism_mapper:  # current precursor ends here
                    assert len(mature_accessions) == len(mature_identifiers)
                    assert precursor_accession
                    assert precursor_identifier

                    taxonomy_id = organism_mapper[organism]
                    mir_precursor_to_mature_file.write('%s\t%s\t%s\n' %
                                                       (taxonomy_id, precursor_identifier,
                                                        ';'.join(mature_identifiers)))
                    organisms_found.add(organism)
                else:
                    organisms_not_mapped.add(organism)

                precursor_identifier = None
                precursor_accession = None
                mature_accessions = []
                mature_identifiers = []
    logger.info('Following organism could not be mapped to taxonomy ID: %s\n' % ','.join(sorted(organisms_not_mapped)))
    logger.debug('\n'.join(organisms_found))


def __get_mir_precursor_to_mature_mapper__():
    """
    :return: a dictionary mapping taxonomy ID and miRBase precursor ID to a list containing the IDs of mature miRNAs
             made from this precursor, e.g. the tuple ('9606', 'hsa-mir-155') might be mapped to ['hsa-miR-155-5p',
             'hsa-miR-155-3p']
    """
    mir_precursor_to_mature = {}
    if logger.getEffectiveLevel() == logging.DEBUG or not os.path.exists(OUT_MIR_PRECURSOR_TO_MATURE_FILE):
        __create_mir_precursor_to_mature_mapping_file__()
    with gzip.open(OUT_MIR_PRECURSOR_TO_MATURE_FILE) as mir_precursor_to_mature_file:
        for line in mir_precursor_to_mature_file:
            taxonomy_id, precursor_identifier, mature_identifiers = line.rstrip('\r\n').split('\t')
            mature_identifiers_split = mature_identifiers.split(';')
            mir_precursor_to_mature[(taxonomy_id, precursor_identifier)] = mature_identifiers_split
    logger.info('Found %i miRNA precursors.' % len(mir_precursor_to_mature))
    for precursor_information, mature_mirs in mir_precursor_to_mature.items():
        taxon_id, mir_precursor = precursor_information
        logger.debug('%s --> %s' % (mir_precursor, ';'.join(mature_mirs)))
    return mir_precursor_to_mature


def create_mir_ortholog_lists():
    """
    Determines groups of miRNA ortholog by following the following steps:
    1) miRBase family annotation (miFam.dat.gz) determines which miRNA precursors belong to the one orthologous group
    2) miRNA precursors in each group are replaced by the mature miRNAs processed from them
    3) only those mature miRNAs that are contained in set of RAIN identifiers are retained
    4) 5prime and 3prime mature miRNAs in the same orthologous group are sorted into two classes

    :return: a dict mapping a string to a list of sets - keys in the returned dictionary are names of miRNA orthologous
             groups and value are a list of sets were set contains the miRBase identifiers of mature miRNAs that are
             orthologous
    """
    # used to check condition 2)
    unique_mir_mapper = stringrnautils.get_unique_mir_mapper()

    # step 1)
    # a list of miRNA families where each family is represented by a set of miRNA precursor names, e.g.,
    # [{'hsa-mir-17', 'hsa-mir-18a', ...} , {'cel-let-7', 'hsa-let-7a-1', 'hsa-let-7a-2' , ...}, .., ]
    precursor_mir_families_list = []
    # a list of miRNA family accessions
    mir_family_accessions = []
    mirs_already_added_to_family = set()  # for sanity checking if a mir precursor is assigned to more than one family
    with gzip.open(MIR_FAMILY_FILE, 'rb') as mir_fam_file:
        current_family_members = set()
        fam_accession = None
        for line in mir_fam_file:
            cols = line.rstrip('\r\n').split()
            tag = cols[0]
            if tag == 'AC':  # family accession
                assert len(cols) == 2
                fam_accession = cols[1]
            elif tag == 'ID':  # family identifier
                assert len(cols) == 2
                # fam_identifier = cols[1]
            elif tag == 'MI':  # accession and identifier of family member
                assert len(cols) == 3
                # fam_member_accession = cols[1]
                fam_member_identifier = cols[2]
                current_family_members.add(fam_member_identifier)
            elif line.startswith('//'):  # end of family
                # sanity check
                for mir in current_family_members:
                    if mir in mirs_already_added_to_family:
                        raise ValueError('miRNA precursor %s is assigned to more than one miRBase family.' % mir)
                mirs_already_added_to_family.update(current_family_members)
                precursor_mir_families_list.append(current_family_members)
                if fam_accession:
                    mir_family_accessions.append(fam_accession)
                else:
                    raise ValueError('Missing family accession.')
                current_family_members = set()
                fam_accession = None
            else:
                raise ValueError('Could not identify tag %s in line %s of file %s.' % (tag, line, MIR_FAMILY_FILE))
    logger.info('Found %i miRNA families in miRBase %s containing a total of %i miRNAs.' %
                (len(precursor_mir_families_list), MIRBASE_VERSION, len(mirs_already_added_to_family)))

    # maps each taxonomy ID and precursor identifier to the identifiers of the mature miRNAs being made from it, e.g.,:
    # mir_precursor_to_mature_mapper[('9606', 'hsa-mir-17')] -> ['hsa-miR-17-5p', 'hsa-miR-17-3p']
    taxonomy_mir_precursor_to_mature_mapper = __get_mir_precursor_to_mature_mapper__()

    # maps each mature miRNA to its taxonomy ID and is needed to produce the final species aware ortholog list. Example:
    # taxonomy_mir_precursor_to_mature_mapper['hsa-miR-17-5p'] -> '9606'
    mature_to_taxonomy_mapper = {}
    # precursors in miR families are not linked to taxonomy ID, so create mapping dictionary of following form:
    # mir_precursor_to_mature_mapper['hsa-mir-17'] -> ['hsa-miR-17-5p', 'hsa-miR-17-3p']
    mir_precursor_to_mature_mapper = {}
    for taxonomy_precursor, mature_list in taxonomy_mir_precursor_to_mature_mapper.items():
        taxonomy_id, precursor_id = taxonomy_precursor
        for mature_id in mature_list:
            mature_to_taxonomy_mapper[mature_id] = taxonomy_id
        mir_precursor_to_mature_mapper[precursor_id] = mature_list

    # step 2)
    mir_precursors_not_mapped = set()
    mature_mir_families_list = []
    for mir_family in precursor_mir_families_list:
        mature_mirs_in_family = set()
        for mir_precursor in mir_family:
            # step 3)
            if mir_precursor in mir_precursor_to_mature_mapper:
                mature_mirs = mir_precursor_to_mature_mapper[mir_precursor]
                for mature_mir in mature_mirs:
                    mature_mirs_in_family.add(mature_mir)
            else:
                mir_precursors_not_mapped.add(mir_precursor)
        mature_mir_families_list.append(mature_mirs_in_family)
    logger.debug('Following miRNAs could not be mapped to mature miRNA(s): ' +
                 ', '.join(sorted(list(mir_precursors_not_mapped))))

    # step 4)
    ends_in_5p = set()
    ends_in_3p = set()
    ends_not_in_3p_or_5p = set()
    final_mir_identifiers = unique_mir_mapper.values()

    # replace each mature miRNA ID by a tuple containing its taxonomy ID and the mature miRNA ID itself
    taxonomy_five_prime_mature_mir_families_list = []
    taxonomy_three_prime_mature_mir_families_list = []
    fam_members_not_mapped = set()
    for family in mature_mir_families_list:
        five_prime_fam_members = set()
        three_prime_fam_members = set()
        for fam_member in family:
            if fam_member not in mature_to_taxonomy_mapper:
                raise IOError('Mature miRNA could not be mapped to taxonomy ID: ' + fam_member)
            taxonomy_id = mature_to_taxonomy_mapper[fam_member]

            if fam_member not in final_mir_identifiers:
                fam_members_not_mapped.add(fam_member)
            elif fam_member.endswith('-5p'):
                five_prime_fam_members.add((taxonomy_id, fam_member))
                ends_in_5p.add(fam_member)
            elif fam_member.endswith('-3p'):
                three_prime_fam_members.add((taxonomy_id, fam_member))
                ends_in_3p.add(fam_member)
            else:
                ends_not_in_3p_or_5p.add(fam_member)
        taxonomy_five_prime_mature_mir_families_list.append(five_prime_fam_members)
        taxonomy_three_prime_mature_mir_families_list.append(three_prime_fam_members)
    logger.info('Following miRNAs were not found in RAIN identifiers: ' + ', '.join(sorted(fam_members_not_mapped)))
    logger.info('Mature miRNAs ending in -5p: %i' % len(ends_in_5p))
    logger.info('Mature miRNAs ending in -3p: %i' % len(ends_in_3p))
    logger.debug('Following miRNAs did neither end in -5p nor on -3p: ' + ', '.join(ends_not_in_3p_or_5p))

    assert (2 * len(mir_family_accessions)) == \
           (len(taxonomy_five_prime_mature_mir_families_list) + len(taxonomy_three_prime_mature_mir_families_list))
    assert len(taxonomy_five_prime_mature_mir_families_list) == len(taxonomy_three_prime_mature_mir_families_list)

    # finally map all family accessions (split into 5p and 3p) families to the mature miRs belonging to this family.
    # Example: fam_map['MIPF0000001-5p'] = set(('9606', 'hsa-miR-17-5p'), ... )
    fam_map = {}
    for i, mir_fam_acc in enumerate(mir_family_accessions):
        fam_map[mir_fam_acc + '-5p'] = taxonomy_five_prime_mature_mir_families_list[i]
        fam_map[mir_fam_acc + '-3p'] = taxonomy_three_prime_mature_mir_families_list[i]
    return fam_map
