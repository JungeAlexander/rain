import gzip
import argparse
import collections
import os
import prepare_miRNA_text_mining_dictionaries as miRNAdictionaries
import logging

logger = logging.getLogger(os.path.basename(__file__))
logging.basicConfig(level=logging.INFO)

parser = argparse.ArgumentParser(prog='PROG',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''
     Generate ortholog file, name file and entity file
     for non-coding (non-miR) RNA.
     These are to be used by LJJ text-mining algorithm.

     Authors: Christian Garde, garde@cbs.dtu.dk, CBS, DTU Systems Biology
              Alexander Junge, RTH, University of Copenhagen <alexander.junge@gmail.com>
                                 ''')

parser.add_argument('-start_serial',
                    type=int,
                    default=100000000)  # LJJ said to use this number instead :D

args = parser.parse_args()

# -----------------------------
# Sub-routines
#-----------------------------


def generate_namefile(alias2name, namefile, start_child_serial):
    """
    Writes text mining name files, a tab-delimited text file with two columns:
    <serial number>\t<name>
    :param alias2name: a dictionary of the form (taxonomy ID, alias) -> identifier
    :param namefile: the path of the name file to be written
    :param start_child_serial: the first serial that is assigned
    :return: a tuple of dictionary and integer; the dictionary is of the form (taxonomy ID, identifier) -> serial number
                                                The integer is the first free serial that can be assigned subsequently.
    """
    # use input to generate dictionary of the form (taxonomy ID, identifier) -> aliases (for identifier)
    name2aliases = collections.defaultdict(set)
    for alias, name in alias2name.iteritems():
        # do not add self-mappings for identifiers here
        if not alias[1] == name:
            name2aliases[(alias[0], name)].add(alias[1])

    # maps tuples of taxon ID and identifier to its integer serial
    name2serial = {}

    with open(namefile, 'w') as handle:
        serial = start_child_serial
        # assign serial to all real names
        for name, aliases in name2aliases.iteritems():
            name2serial[name] = serial
            handle.write("%d\t%s\n" % (serial, name[1]))
            for alias in aliases:
                handle.write("%d\t%s\n" % (serial, alias))
            serial += 1

    return name2serial, serial



def generate_groupfile(orth_dict, name2serial, groupfile, start_parent_serial):

    with open(groupfile, 'w') as handle:
        serial = start_parent_serial
        group2serial = {}
        for key, vals in orth_dict.iteritems():
            handle.write("%d\t%d\n" % (name2serial[key], serial))
            group2serial["OrthGroup_%s" % ":".join(key)] = serial
            for val in vals:
                handle.write("%d\t%d\n" % (name2serial[val], serial))
            serial += 1

        all_orth = set([z for x, y in orth_dict.iteritems() for z in [x] + y])
        for key in name2serial.keys():
            if not key in all_orth:
                handle.write("%d\t%d\n" % (name2serial[key], serial))
                group2serial["OrthGroup_%s" % ":".join(key)] = serial
                serial += 1

    return serial, group2serial


def generate_mirna_groupfile(mirna_ortholog_dict, name2serial, groupfile, start_parent_serial):
    """
    Writes group file storing information about orthologous miRNAs that us uses during text mining.
    The group file is tab-delimited and hast two columns:
    <serial of member of ortholog group>\t<serial of ortholog group identifier>

    :param mirna_ortholog_dict: Maps name of a miRNA ortholog group to a set containing tuples of taxonomy ID and name
                                of mature miRNAs in that group. E.g.:
                                mir_text_mining_ortholog_groups['MIPF0001322-3p'] =
                                                                {('7955', 'dre-miR-459-3p'), ('8090', 'ola-miR-459-3p')}
    :param name2serial: dictionary of the form (taxonomy ID, miR identifier) -> serial number
    :param groupfile: path of the group file to be written
    :param start_parent_serial: the first serial that is assigned
    :return: tuple of integer and dictionary - Integer is first free serial that can be assigned in subsequent steps,
                                               dictionary maps group identifiers to serial, e.g., ['MIPF0001322-3p']->12
    """
    serial = start_parent_serial
    group2serial = {}
    identifiers_in_ortho_group = set()  # set of tuples (tax, miR identifier), all contained tuples are in orthol. group
    with open(groupfile, 'w') as handle:
        for mir_text_mining_group_id, mir_text_mining_group_members_set in mirna_ortholog_dict.iteritems():
            # Skip empty miRNA ortholog groups
            if len(mir_text_mining_group_members_set) == 0:
                continue
            # assign serial to current group
            group2serial[mir_text_mining_group_id] = serial
            # associate all orthologs' serials with group serial in ortholog file
            for tax_miR_identifier in mir_text_mining_group_members_set:
                if tax_miR_identifier in identifiers_in_ortho_group:
                    raise IOError('Following tuple seems to be in one than more ortholog group.' +
                                  str(tax_miR_identifier))
                identifiers_in_ortho_group.add(tax_miR_identifier)
                handle.write('%d\t%d\n' % (name2serial[tax_miR_identifier], serial))
            serial += 1

        # Create singleton orthologous group for each identifier that has not been assigned to any group so far
        for key in name2serial.keys():
            if not key in identifiers_in_ortho_group:
                handle.write("%d\t%d\n" % (name2serial[key], serial))
                group2serial["OrthGroup_%s" % ":".join(key)] = serial
                serial += 1
    return serial, group2serial


def generate_entityfile(name2serial, group2serial, entityfile):
    with open(entityfile, 'w') as handle:
        for tax_identifier, identifier_serial in name2serial.iteritems():
            tax, identifier = tax_identifier
            handle.write("%d\t%s\t%s\n" % (identifier_serial, tax, identifier))
        for group_id, group_serial in group2serial.iteritems():
            # hardcoded! 2759 is the type of our groups, because all our organisms are eukaryotes
            # ideally we would recode this to be part of the dictionary so we could set it to 1 if we
            # include bacteria or archaea later
            group_type = "2759"
            handle.write("%d\t%s\t%s\n" % (group_serial, group_type, group_id))


def rm_dubious_alias(orth_dict):
    # Remove aliases that point to multiple prime ids
    aliases = [[i, y] for i, (x1, x2) in enumerate(orth_dict.iteritems()) for y in x2]
    alias_lookup = collections.defaultdict(list)
    for idx, alias in aliases:
        alias_lookup[alias].append(idx)
    alias2rm = set([alias for alias, idxes in alias_lookup.iteritems() if len(idxes) > 1])
    new_orth = {}
    for name, alias in orth_dict.iteritems():
        new_orth[name] = [alias for alias in orth_dict[name] if not alias in alias2rm]
    return new_orth


def ncrna_generate_orthlist(orth_file, alias2name):
    ## Exclude database identifiers? Not included currently
    #exclude = set([ "ENSG", "ENSM", "ENSDARG", "NR_", "HGNC:",
    #       "FBan", "FBgn", "WBGene","ZDB-GENE-","ENSOCUG","ENST",
    #       "MGI:"])

    # Taxonomy in orth file: human, mouse, zebrafish, fly, worm, hard coding :(
    taxon = ['9606', '9606', '10090', '10090', '7955', '7955', '7227', '7227', '6239', '6239']

    # Load and filter orthologues to those included in ncRNA alias file (e.g. exclude miR and proteins)
    handle = gzip.open(orth_file) if orth_file.endswith(".gz") else open(orth_file)
    orth = [x.strip("\n").split("\t") for x in handle]
    handle.close()
    orth = [x for x in orth if ('9606', x[1]) in alias2name]

    # Generate orth registry
    orth_dict = collections.defaultdict(set)
    for record in orth:
        if len(record[1]) == 0:
            continue
        key = ("9606", record[1])
        for idx in range(2, 10, 2):  # some hard coding :(
            if record[idx + 1] == "1" and (taxon[idx], record[idx]) in alias2name:
                orth_dict[(key[0],alias2name[key])].add((taxon[idx], alias2name[(taxon[idx], record[idx])]))

    orth_dict = rm_dubious_alias(orth_dict)
    return orth_dict


def retrieve_ncrna_files(alias_file):
    # Load alias convertion, format: tax (tab) official name (tab) alias (tab) database
    handle = gzip.open(alias_file) if alias_file.endswith(".gz") else open(alias_file)
    alias2names = [x.strip().split("\t") for x in handle]
    handle.close()
    alias2name = dict([((x[0], x[2]), x[1]) for x in alias2names if not x[3] == "Ensembl_Description"]) # remove the description, as that probably don't work well for the tagger
    return alias2name


def main():
    odir = "id_dictionaries"

    #----------------
    # ncRNA
    #----------------
    ncrna_namefile = os.path.join(odir, "ncRNA_namefile.tsv")
    ncrna_entityfile = os.path.join(odir, "ncRNA_entityfile.tsv")
    ncrna_groupfile = os.path.join(odir, "ncRNA_groupfile.tsv")
    ncrna_aliasfile = os.path.join(odir, "ncRNAaliasfile.tsv.gz")
    ncrna_orthfile = os.path.join(odir, "ncRNAorthfile.tsv.gz")

    alias2name = retrieve_ncrna_files(ncrna_aliasfile)

    # Generate textmining files
    ncrna_orthdict = ncrna_generate_orthlist(ncrna_orthfile, alias2name)
    name2serial, start_parent_serial = generate_namefile(alias2name, ncrna_namefile, args.start_serial)
    start_mir_serial, group2serial = generate_groupfile(ncrna_orthdict, name2serial, ncrna_groupfile,
                                                        start_parent_serial)
    generate_entityfile(name2serial, group2serial, ncrna_entityfile)

    # start_mir_serial = args.start_serial # debugging tingy :)
    #----------------
    # miRNA - implemented by Alexander Junge, RTH, University of Copenhagen <alexander.junge@gmail.com>
    #----------------
    mir_namefile = os.path.join(odir, 'miR_namefile.tsv')
    mir_entityfile = os.path.join(odir, 'miR_entityfile.tsv')
    mir_groupfile = os.path.join(odir, 'miR_groupfile.tsv')

    # Maps a tuple of taxonomy ID and miRNA alias to final miRNA identifier.
    # Example: mir_text_mining_alias_dictionary[('7739', 'MIMAT0020086')] = 'bfl-miR-4859-5p'
    mir_text_mining_alias_dictionary = miRNAdictionaries.get_text_mining_mir_dictionary()

    # Maps name of a miRNA ortholog group to a set containing tuples of taxon ID and name of mature miRNAs in that group
    # E.g.: mir_text_mining_ortholog_groups['MIPF0001322-3p'] =,{('7955', 'dre-miR-459-3p'), ('8090', 'ola-miR-459-3p')}
    mir_text_mining_ortholog_groups = miRNAdictionaries.create_mir_ortholog_lists()

    name2serial, start_mirparent_serial = generate_namefile(mir_text_mining_alias_dictionary, mir_namefile,
                                                            start_mir_serial)
    trash_no, group2serial = generate_mirna_groupfile(mir_text_mining_ortholog_groups, name2serial, mir_groupfile,
                                                      start_mirparent_serial)
    generate_entityfile(name2serial, group2serial, mir_entityfile)


if __name__ == '__main__':
    logger.info('Generating ortholog file, name file and entity file files needed for text mining.')
    main()
    logger.info('Done.' + os.linesep + os.linesep)