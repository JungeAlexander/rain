import stringrnautils
import os, argparse, gzip
import logging
logger = logging.getLogger(os.path.basename(__file__))
logging.basicConfig(level=logging.INFO)

def is_valid_path(arg):
    if not os.path.exists(arg):
        parser.error("Given master file path %s does not exist." % arg)
    return arg


parser = argparse.ArgumentParser(description=
                                 """Process Lars Juhl Jensen's Textmining.
                                    Author: Christian Garde, Workman Lab, DTU Systems Biology;
                                            Alexander Junge, RTH, University of Copenhagen""")
parser.add_argument('gold_standard_file', type=is_valid_path,
                    help='The gold standard file to benchmark against.')
parser.add_argument('-master_path',
                    help='path to write master files',
                    default = 'master_files')
parser.add_argument('-data_path',
                    help='path to write data files',
                    default = 'data')
parser.add_argument('-id_path',
                    help='path to write id dictionaries',
                    default = 'id_dictionaries')
parser.add_argument('-disallow_mirpairs',
                    help='whether to disallow miR pairs',
                    action='store_true')
parser.add_argument('-disallow_sno_mir_pairs',
                    help='whether to disallow snoRNA-miR pairs',
                    action='store_true')
parser.add_argument('-verbose',
                    help='whether to report when identifiers cannot be mapped',
                    action='store_true')
logger.info('Parsing Textmining')
args = parser.parse_args()
gold_standard_file_path = args.gold_standard_file
id_dict_dir = args.id_path
data_dir = args.data_path
master_dir = args.master_path
disallow_mirpairs = args.disallow_mirpairs
disallow_sno_mir_pairs = args.disallow_sno_mir_pairs

if disallow_mirpairs:
    logger.info('Removing miRNA-miRNA interactions.')

if disallow_sno_mir_pairs:
    logger.info('Removing snoRNA-miRNA interactions.')


def open_file(x):
    return gzip.open(x) if x.endswith(".gz") else open(x)

# textmining results
textmining_file = os.path.join(data_dir, "rain_textmining.tsv.gz")

# load ncRNA dictionary
ncrna_mapper = stringrnautils.get_non_coding_rna_alias_mapper()

# Load miR mapper
mir_mapper = stringrnautils.get_alias_mir_mapper()

# snoRNA dictionary
snoRNA_file = os.path.join(id_dict_dir, "snoRNA_lookup.tsv.gz")

snoRNA = [ x.strip("\n").split("\t") for x in open_file(snoRNA_file) ]
snoRNA = [ (x[0], ncrna_mapper[x[0]][x[1]] ) for x in snoRNA if x[1] in ncrna_mapper[x[0]]]

# Load string v.9 --> string v10 mapper...supposedly the text mining was run with string v9 identifiers
# quick fix for human, other organisms didn't have the problem...
ENSP2ENSG_v8  = stringrnautils.get_string_to_alias_mapper('9606','ENSP','ENSG', 8, 'all', True)['9606']
ENSP2ENSG_v9  = stringrnautils.get_string_to_alias_mapper('9606','ENSP','ENSG', 9, 'all', True)['9606']
ENSG2ENSP_v10  = stringrnautils.get_alias_to_string_mapper(9606, 'ENSP', 'ENSG',10 )['9606']
ENSP8_to_ENSP10 = dict([ (ensp,ENSG2ENSP_v10[ensg]) for ensp,ensg in ENSP2ENSG_v8.iteritems() if ensg in ENSG2ENSP_v10 ])
ENSP9_to_ENSP10 = dict([ (ensp,ENSG2ENSP_v10[ensg]) for ensp,ensg in ENSP2ENSG_v9.iteritems() if ensg in ENSG2ENSP_v10 ])
ENSP2ENSP_v10 = stringrnautils.get_alias_to_string_mapper(9606, 'ENSP', 'ENSP',10 )['9606']

old_ensp2ncrna = dict([(ensp, ncrna_mapper['9606'][ensg]) for ensp,ensg in ENSP2ENSG_v8.iteritems()
                       if ensg in ncrna_mapper['9606']])

old_ensp2ncrna.update(dict([(ensp, ncrna_mapper['9606'][ensg]) for ensp, ensg in ENSP2ENSG_v9.iteritems()
                       if ensg in ncrna_mapper['9606']]))

# Benchmark textmining
def check_for_ensp_issue(ID, ENSP2ENSG_v8, ENSP2ENSG_v9, ENSP8_to_ENSP10, ENSP9_to_ENSP10, ENSP2ENSP_v10,
                         old_ensp2ncrna):
    is_protein = False
    if ID in ENSP2ENSP_v10:
        ID = ENSP2ENSP_v10[ID]
        is_protein = True
    elif ID in ENSP2ENSG_v9:
        if ID in ENSP9_to_ENSP10:
            ID = ENSP9_to_ENSP10[ID]
            is_protein = True
        elif ID in old_ensp2ncrna:
            ID = old_ensp2ncrna[ID]
            is_protein = True
        else:
            logger.debug("Unable to convert %s from v9 to v10 of string" % ID)
            ID = -1
    elif ID in ENSP2ENSG_v8:
        if ID in ENSP8_to_ENSP10:
            ID = ENSP8_to_ENSP10[ID]
            is_protein = True
        elif ID in old_ensp2ncrna:
            ID = old_ensp2ncrna[ID]
            is_protein = True
        else:
            logger.debug("Unable to convert %s from v8 to v10 of string" % ID)
            ID = -1
    elif ID.startswith("ENSP"):
        logger.debug("%s neither found in v8, v9 or v10 of string" % ID)
        ID = -1
    else:
        pass
    return is_protein, ID

def check_sno_miR(record,mir_mapper, snoRNA):
    if record[1] in mir_mapper and (record[0],record[3]) in snoRNA:
        is_pair=True
    elif record[3] in mir_mapper and (record[0],record[1]) in snoRNA:
        is_pair=True
    else:
        is_pair=False
    return is_pair

def check_miR_pair(record,mir_mapper):
    return True if record[1] in mir_mapper and record[3] in mir_mapper else False

logger.info('Reading text mining interactions.')
pairs = {}
string_species = stringrnautils.get_string_10_species()
protein_protein_pairs = set()
for record in gzip.open(textmining_file):
    # format (tab-delimited):
    # taxonomy_id   ID1   alias_human_readable1  ID2 alias_human_readable2  z-score stars   link
    score_idx = 5  # z-scores
    record = record.strip("\n").split("\t")

    if record[0] not in string_species:  # reduce to organisms of interest
        continue

    first_is_protein, record[1] = check_for_ensp_issue( record[1], ENSP2ENSG_v8, ENSP2ENSG_v9,
                                      ENSP8_to_ENSP10, ENSP9_to_ENSP10, ENSP2ENSP_v10, old_ensp2ncrna)
    if record[1] == -1:
        continue

    second_is_protein, record[3] = check_for_ensp_issue( record[3], ENSP2ENSG_v8, ENSP2ENSG_v9,
                                      ENSP8_to_ENSP10, ENSP9_to_ENSP10, ENSP2ENSP_v10, old_ensp2ncrna)

    if record[3] == -1:
        continue

    key = ( record[0], record[1], record[3] )
    rev_key = ( record[0], record[3], record[1] )

    if first_is_protein and second_is_protein:
        if key in protein_protein_pairs or rev_key in protein_protein_pairs:
            pass
        else:
            protein_protein_pairs.add(key)
        continue

    if key in pairs or rev_key in pairs: # don't include duplicates
        continue
    elif disallow_mirpairs and check_miR_pair(record,mir_mapper): # exclude miR pairs
        continue
    elif disallow_sno_mir_pairs and check_sno_miR(record,mir_mapper,snoRNA): # exclude sno-miR pairs
        continue
    elif record[3] in mir_mapper: # always have miR in first position for benchmarking
        pairs[rev_key] = float( record[score_idx] )
    elif record[1] in ncrna_mapper[record[0]] or record[1] in mir_mapper:
        pairs[key] = float(record[score_idx])
    else:
        pass

logger.info('Removed {:d} protein-protein interactions from RAIN text mining results'.format(
    len(protein_protein_pairs)))
orgn,rna,other,scores = zip(*[ (x[0],x[1],x[2],y) for x,y in pairs.iteritems() ])

logger.info('Benchmarking.')
new_scores = stringrnautils.benchmark(orgn, rna, other, scores, gold_standard_file_path,
                                          increases=True, window_size=100,
                                          fit_name="textmining")

# Write to master file
data_out = [ "\t".join((tax,entity1,entity2,"0","textmining",str(score), "Textmining", "", ""))+'\n'
             for tax,entity1,entity2,score in zip(orgn,rna,other,new_scores) ]

logger.info('Writing master file.')
with open(os.path.join(master_dir, "textmining.tsv"),'w') as handle:
    handle.write("".join(data_out))

logger.info('Done.' + os.linesep + os.linesep)

