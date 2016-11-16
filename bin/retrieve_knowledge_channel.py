import os
import argparse
import tempfile
import shutil
import re
import stringrnautils
import logging

logger = logging.getLogger(os.path.basename(__file__))
logging.basicConfig(level=logging.INFO)

parser = argparse.ArgumentParser(description=
                                 """Retrieves manually curated knowledge channel.

                                    Authors: Christian Garde, Workman Lab, DTU Systems Biology
                                             Alexander Junge, RTH, University of Copenhagen""")
parser.add_argument('-master_path',
                    help='path to write master files',
                    default='master_files')
parser.add_argument('-data_path',
                    help='path to write data files',
                    default='data')
args = parser.parse_args()
DATA_PATH = args.data_path
MASTER_PATH = args.master_path

logger.info('Retrieving Knowledge channel. Also creating data and master file directories if they do not exist.')

# Final version
knowledge_file = os.path.join(DATA_PATH, 'knowledge_channel.tar.gz')

if not os.path.exists(DATA_PATH):
    os.mkdir(DATA_PATH)

if not os.path.exists(MASTER_PATH):
    os.mkdir(MASTER_PATH)

cat_master_file_name = 'database.tsv'
knowledge_master_file_path = os.path.join(MASTER_PATH, cat_master_file_name)

################################################################################
# NOTE/TODO for future versions of RAIN:
# currently, this script assumes that the curated knowledge channel contains
# only human interactions. If this changes in the future, the ID mapping
# routines must be adjusted.
################################################################################

# Make up for the fact that it was curated for version 9
ENSP2ENSG_v9  = stringrnautils.get_string_to_alias_mapper('9606','ENSP','ENSG', 9, 'all', True)['9606']
ENSG2ENSP_v10  = stringrnautils.get_alias_to_string_mapper(9606, 'ENSP', 'ENSG',10 )['9606']
ENSP9_to_ENSP10 = dict([(ensp,ENSG2ENSP_v10[ensg]) for ensp,ensg in ENSP2ENSG_v9.iteritems() if ensg in ENSG2ENSP_v10])
ENSP9_to_ENSP10.update({'ENSP00000403359': 'ENSP00000441000', 'ENSP00000400867': 'ENSP00000441000',
                        'ENSP00000403175': 'ENSP00000393241'})

ENSP2ENSP_v10 = stringrnautils.get_alias_to_string_mapper(9606, 'ENSP', 'ENSP',10 )['9606']
ncrna_mapper = stringrnautils.get_non_coding_rna_alias_mapper()['9606']

def correct_rna_names(ID, ncrna_mapper):
    if ID == 'RNA5-8S1':
        ID = "RNA5-8S5"
    elif ID == 'RNA28S1':
        ID = 'RNA28S5'
    elif ID == 'RNA18S1':
        ID = 'RNA45S5'
    if ID in ncrna_mapper:
        ID = ncrna_mapper[ID]
    return ID

def check_for_ensp_issue( ID, ENSP2ENSG_v9, ENSP9_to_ENSP10, ENSP2ENSP_v10, verbose=True ):
    failed_pairs = []
    if ID in ENSP2ENSP_v10:
        ID = ENSP2ENSP_v10[ID]
    elif ID in ENSP2ENSG_v9:
        if ID in ENSP9_to_ENSP10:
            ID = ENSP9_to_ENSP10[ID]
        else:
            logger.warning('Unable to convert %s (ensg: %s ) from v9 to v10 of STRING.' % (ID, ENSP2ENSG_v9[ID]))
            return False
    return ID

output=[]
if not os.path.exists(knowledge_master_file_path):
    # extract master files for single RNA classes to temporary directory
    tmp_dir_path = tempfile.mkdtemp()
    os.system("tar xzf %s -C %s" % (knowledge_file, tmp_dir_path))

    # cat all master files for single RNA classes to one master file
    single_rna_class_master_files = os.listdir(tmp_dir_path)
    cat_master_file_path = os.path.join(tmp_dir_path, cat_master_file_name)
    
    # remove potentially duplicated interaction lines
    interactions_seen = set()
    for single_class_file_path in single_rna_class_master_files:
        logger.info('Adding knowledge file %s.' % single_class_file_path)
        with open(os.path.join(tmp_dir_path, single_class_file_path), 'r') as single_class_file:
            for line in single_class_file:
                if re.match('^\s*#', line):
                    continue
                row = line.rstrip('\r\n').split('\t')
                try:
                    organism, node1, node2, directed, evidence, score, source, linkout, comment = row
                    if organism != '9606':
                        continue
                    node1 = correct_rna_names( node1,ncrna_mapper)
                    node2 = correct_rna_names( node2,ncrna_mapper)
                    node1 = check_for_ensp_issue( node1, ENSP2ENSG_v9, ENSP9_to_ENSP10,ENSP2ENSP_v10)
                    node2 = check_for_ensp_issue( node2, ENSP2ENSG_v9, ENSP9_to_ENSP10,ENSP2ENSP_v10)
                    if node1==False or node2==False:
                        continue
                except ValueError as e:
                    logger.critical('Error (see below) encountered when unpacking row in aggregate master file. '
                                    'Corresponding row was:\n%s\n' % row)
                    raise e
                interaction_key = (organism, node1, node2)
                if interaction_key in interactions_seen:
                    continue
                interactions_seen.add(interaction_key)
                # fix evidence  to 'database'
                evidence = 'database'
                output.append('\t'.join((organism, node1, node2, directed, evidence, score, source,
                                                 linkout, comment))+'\n' )

    with open(cat_master_file_path, 'w') as cat_master_file:
        cat_master_file.write("".join([ x for x in set(output) ]) )

    # copy concatenated file to master file
    shutil.copy(cat_master_file_path, MASTER_PATH)

    # remove temporary directory
    shutil.rmtree(tmp_dir_path)
logger.info('Done.' + os.linesep + os.linesep)
