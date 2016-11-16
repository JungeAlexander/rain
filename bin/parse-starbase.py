import os, argparse, gzip
import stringrnautils
import logging

logger = logging.getLogger(os.path.basename(__file__))
logging.basicConfig(level=logging.INFO)

def is_valid_path(arg):
    if not os.path.exists(arg):
        parser.error("Given master file path %s does not exist." % arg)
    return arg


parser = argparse.ArgumentParser(description=
                                 """Process starbase, clash and crac studies.
                                    Author: Christian Garde, Workman Lab, DTU Systems Biology""")
parser.add_argument('gold_standard_file', type=is_valid_path,
                    help='The gold standard file to benchmark against.')
parser.add_argument('-master_path',
                    help='path to write master files',
                    default = 'master_files')
parser.add_argument('-data_path',
                    help='path to write data files',
                    default = 'data')
parser.add_argument('-id_path',
                    help='path to write data files',
                    default = 'id_dictionaries')

args = parser.parse_args()

logger.info('Parsing starbase, clash and crac.')

logger.debug(
    """

    PROCESSING STARBASE 2.0, CLASH and CRAC STUDIES
    ------------------------------------------------
    Last revision: Nov 3, 2014

    Resources:
    [1] STARBASE 2.0: NAR 2014;42:D92-7. doi: 10.1093
    [2] HUMAN CLASH: Cell 2013;153(3):654-65. doi: 10.1016
    [3] YEAST CLASH: PNAS 2011;108(24):10010-5. doi: 10.1073
    [4] YEAST CRAC: Cell 2013;154(5):996-1009. doi: 10.1016

    Christian Garde
    """)

# "https://docs.google.com/uc?authuser\u003d0\u0026id\u003d0ByN5rGHG0ELyZi16a05OdHBwU0E\u0026export\u003ddownload".decode("unicode_escape")

"""
Reference URLS:
HUMAN CLASH: http://www.cell.com/cms/attachment/2007954362/2030562248/mmc1.txt / http://www.ncbi.nlm.nih.gov/pubmed/23622248
YEAST CLASH: http://www.pnas.org/content/suppl/2011/05/24/1017386108.DCSupplemental/sd01.xls
STARBASE: http://starbase.sysu.edu.cn/
CRAC: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46742 / http://www.ncbi.nlm.nih.gov/pubmed/23993093

Fetching data from public google drive copy to avoid risk of respective servers being unavailable
"""

##########################
# Generate Directory Structure
##########################

id_dict_dir = args.id_path
data_dir = args.data_path
master_dir = args.master_path
gold_standard_file_path = args.gold_standard_file

print_unknown = 0

##########################
# Extract data data
##########################


starbase_file = os.path.join(data_dir, "starbase.tar.gz")
clash_file = os.path.join(data_dir,"clash.tar.gz")
crac_file = os.path.join(data_dir, "crac.tsv.gz")
ncrna_file = os.path.join(id_dict_dir,"ncRNAaliasfile.tsv.gz")
mrna_file = os.path.join(id_dict_dir, "ensembl72_protein_lookup.tsv.gz"  )

logger.debug("STEP 1: Extracting starbase data...")
os.system("tar -zxf %s -C %s" % (starbase_file, data_dir))

logger.debug("STEP 2: Retrieving clash data...")
os.system("tar -zxf %s -C %s" % (clash_file, data_dir))

#############################
# Load identifier dictionaries
#############################

logger.debug("STEP 4: Retrieving mapping dictionaries...")

# Load miR mapper
mir_mapper = stringrnautils.get_unique_mir_mapper()

# Load ncRNA mapper
ncrna_mapper = stringrnautils.get_non_coding_rna_alias_mapper()

# Load STRING_mapper
string_mapper = stringrnautils.get_alias_to_string_mapper(9606, 'ENSP', 'ENSG')
string_mapper.update( stringrnautils.get_alias_to_string_mapper('10090', 'ENSMUSP', 'ENSMUSG') )
string_mapper.update( stringrnautils.get_alias_to_string_mapper('6239', '', 'WBGene') )
string_mapper.update( stringrnautils.get_alias_to_string_mapper('4932', '', '') )
string_mapper = dict([ [(taxid,y[0]),y[1]] for taxid,x in string_mapper.iteritems() for y in x.iteritems() ])

handle = gzip.open(mrna_file) if mrna_file.endswith(".gz") else open(mrna_file)
mrna_mapper = [ x.strip("\n").split("\t") for x in handle ]
handle.close()
mrna_mapper = dict([ [(x[0], x[2]),string_mapper[ (x[0], x[1]) ] ] for x in mrna_mapper if (x[0], x[1]) in string_mapper ])

# Load 'hard coded' protein dictionary, because the starbase protein clip names are non-standard
protein_mapper = { ('9606','MOV10')     :'ENSG00000155363',
                   ('9606','FUS')       :'ENSG00000089280',
                   ('9606','eIF4AIII')  :'ENSG00000141543',
                   ('9606','C22ORF28')  :'ENSG00000100220',
                   ('9606','UPF1')      :'ENSG00000005007',
                   ('9606','HuR')       :'ENSG00000066044',
                   ('9606','CAPRIN1')   :'ENSG00000135387',
                   ('9606','C17ORF85')  :'ENSG00000074356',
                   ('9606','DGCR8')     :'ENSG00000128191',
                   ('9606','QKI')       :'ENSG00000112531',
                   ('9606','FUS-mutant'):'ENSG00000089280',
                   ('9606','PTB')       :'ENSG00000011304',
                   ('9606','SFRS1')     :'ENSG00000136450',
                   ('9606','ZC3H7B')    :'ENSG00000100403',
                   ('9606','LIN28')     :'ENSG00000131914',
                   ('9606','LIN28B')    :'ENSG00000187772',
                   ('9606','LIN28A')    :'ENSG00000131914',
                   ('9606','U2AF65')    :'ENSG00000063244',
                   ('9606','FXR1')      :'ENSG00000114416',
                   ('9606','FXR2')      :'ENSG00000129245',
                   ('9606','TDP43')     :'ENSG00000120948',
                   ('9606','IGF2BP1')   :'ENSG00000159217',
                   ('9606','IGF2BP3')   :'ENSG00000136231',
                   ('9606','IGF2BP2')   :'ENSG00000073792',
                   ('9606','TNRC6')     :'ENSG00000090905',
                   ('9606','PUM2')      :'ENSG00000055917',
                   ('9606','TAF15')     :'ENSG00000172660',
                   ('9606','ALKBH5')    :'ENSG00000091542',
                   ('9606','TIA1')      :'ENSG00000116001',
                   ('9606','TIAL1')     :'ENSG00000151923',
                   ('9606','FMRP')      :'ENSG00000102081',
                   ('9606','hnRNPC')    :'ENSG00000092199',
                   ('9606','EWSR1')     :'ENSG00000182944',
                   ('10090','Cirbp')    :'ENSMUSG00000045193',
                   ('10090','EZH2')     :'ENSMUSG00000029687',
                   ('10090','FUS')      :'ENSMUSG00000030795',
                   ('10090','MBNL1')    :'ENSMUSG00000027763',
                   ('10090','MBNL2')    :'ENSMUSG00000022139',
                   ('10090','nElavl')   :'ENSMUSG00000028546',
                   ('10090','NOVA')     :'ENSMUSG00000021047',
                   ('10090','PTBP2')    :'ENSMUSG00000028134',
                   ('10090','Rbm3')     :'ENSMUSG00000031167',
                   ('10090','SRSF1')    :'ENSMUSG00000018379',
                   ('10090','SRSF2')    :'ENSMUSG00000034120',
                   ('10090','TDP43')    :'ENSMUSG00000041459',
                   ('6239','GLD1')      :'WBGene00001595',
                   ('4932','Cbc1')      :'YML020W',
                   ('4932','Gbp2')      :'YCL011C',
                   ('4932','Hek2')      :'YBL032W',
                   ('4932','Hrp1')      :'YOL123W',
                   ('4932','Mex67')     :'YPL169C',
                   ('4932','Mtr4')      :'YJL050W',
                   ('4932','Nab2')      :'YGL122C',
                   ('4932','Pab1')      :'YER165W',
                   ('4932','Ski2')      :'YLR398C',
                   ('4932','Tho2')      :'YNL139C',
                   ('4932','Tif1')      :'YKR059W',
                   ('4932','Trf4')      :'YOL115W',
                   ('4932','Xrn1')      :'YGL173C'
                 }

protein_mapper = dict([ [x, string_mapper[x[0], y]] for x,y in protein_mapper.iteritems() ])

##################################################
# Define Subrutines
###################################################


def convert_id(x,x_type):
    if x_type == 1: # miRNA
        return mir_mapper[x[1]]
    if x_type == 2: # non-miRNA-ncRNA
        return ncrna_mapper[x[0]][x[1]]
    if x_type == 3: # mRNA
        return mrna_mapper[x]
    if x_type == 4: # CLIP protein
        return protein_mapper[x]
    if x_type == 5: # GeneIdentifier --> StringProteinIdentifier
        return string_mapper[x]
    if x_type == -1:
        return x


def get_convert_list( filename, mol1_type, mol2_type, idx1, idx2 ):
    filename_base = filename.split("/")[-1]
    output = []
    tax_id = filename_base.split("_")[0]
    n_exp = int(filename.split("=")[-1].split(".")[0])
    direction = '1' if 'Starbase-miRNA-mRNA' in filename_base else '0'
    slice_name = 1 if 'Starbase-miRNA-circRNA' in filename_base else 0
    no_unmapped = 0
    handle = gzip.open(filename) if filename.endswith(".gz") else open( filename )
    for idx, record in enumerate(handle):
        if idx == 0: # header
            continue
        record = record.split("\t")
        if slice_name:
            record[idx2] = "_".join( record[idx2].split("_")[1:] )
        try:
            id1 = convert_id( (tax_id, record[idx1]), mol1_type)
            id2 = convert_id( (tax_id, record[idx2]), mol2_type)
            output.append([tax_id,id1,id2,direction,'experimental',n_exp,'StarBase2.0', "",""])
        except:
            no_unmapped += 1
            if print_unknown:
                logger.error("Unmapped pair [taxid:%d]: %s, %s\n" % (tax_id, record[idx1],record[idx2]) )
            else:
                pass
    return output


def get_clash_list( filename ):
    filename_base = filename.split("/")[-1]
    output = []
    tax_id = filename_base.split("_")[0]
    direction = '1' if tax_id == '9606' else '0' # human study is miRNA-mRNA whereas yeast study is snoRNA-rRNA
    if tax_id == '9606':
        idx1,idx2,idx_score,mol1_type,mol2_type = (1,5,10,1,5)
        ref_url = 'http://www.ncbi.nlm.nih.gov/pubmed/23622248'
    elif tax_id == '4932':
        idx1,idx2,idx_score,mol1_type,mol2_type = (1,2,6,2,2)
        ref_url = 'http://www.ncbi.nlm.nih.gov/pubmed/21610164'
    else:
        logger.critical("Unknown idx for CLASH study\n")
        exit(1)
    handle = gzip.open(filename) if filename.endswith(".gz") else open( filename )
    for idx, record in enumerate(handle):
        if record.startswith("#") or len(record.strip())==0: # header
            continue
        record = record.strip("\n").split("\t")
        try:
            id1 = convert_id( (tax_id, record[idx1].split("_")[0]), mol1_type)
            id2 = convert_id( (tax_id, record[idx2].split("_")[0]), mol2_type)
            n_exp = int(record[idx_score])
            output.append([tax_id,id1,id2,direction,'experimental', n_exp, 'CLASH', ref_url,""])
        except:
            if print_unknown:
                logger.error("Unmapped pair [taxid:%d]: %s, %s\n" % (tax_id, record[idx1],record[idx2]) )
            else:
                pass
    return output

def get_crac_list( filename ):
    output = []
    tax_id = "4932"
    handle = gzip.open(filename) if filename.endswith(".gz") else open(filename)
    for record in handle:
        record = record.strip("\n").split("\t")
        id1 = convert_id(("4932",record[0]), 4)
        try:
            id2 = convert_id(("4932", record[1]), 2)
        except:
            try:        # try mir-conversion
                id2 = convert_id(record[1], 1)
            except:
                continue # don't try mrna?!?
        output.append( [tax_id, id1, id2, "0", "experimental", int(record[2]), "Crac", "http://www.ncbi.nlm.nih.gov/pubmed/23993093",""] )
    return output

def redundancy_reduce( data ):
    output = {}
    for record in data:
        output_key = tuple(record[:3])
        if output_key in output:
            if record[5] > output[output_key][5]:
                output[output_key] = record
        else:
            output[output_key] = record
    return output.values()

def merge_bin(x,y,larger):
    if larger:
        return str(y) if int(x) > y else str(x)
    else:
        return str(y) if int(x) < y else str(x)

##################################################
# Load Starbase files
###################################################

logger.debug("STEP 5: Converting starbase, clash and crac data to master file format...")

#---------------------
# Co-expression data - determined by hypergeometric test, cf. Starbase v2.0 paper
#---------------------

iFiles = [ os.path.join(data_dir, x) for x in os.listdir(data_dir) if 'Starbase-miRNA' in x and x.endswith(".tsv.gz")]

# Convert ids and format data
coexpr_data = []
idx1_type = 1
for iFile in iFiles:
    idx2_type = 3 if 'Starbase-miRNA-mRNA' in iFile.split("/")[-1] else 2 # define type of molecule 2
    idx1,idx2 = (0,1) if idx2_type==3 else (1,2)
    coexpr_data += get_convert_list( iFile, idx1_type, idx2_type, idx1, idx2 )

coexpr_data = redundancy_reduce(coexpr_data)

#---------------------
# CLIP data
#---------------------

iFiles = [ os.path.join(data_dir, x) for x in os.listdir(data_dir) if 'Starbase-protein' in x and x.endswith(".tsv.gz")]

# Convert ids and format data
clip_data = []
idx1_type = 4
for iFile in iFiles:
    idx2_type = 3 if 'Starbase-protein-mRNA' in iFile.split("/")[-1] else 2 # define type of molecule 2
    if idx2_type == 3: # Skip protein-mRNA interactions? This is mainly splicing events, and according to Lars those should not be included
        continue
    idx1,idx2 = (0,1)
    clip_data += get_convert_list( iFile, idx1_type, idx2_type, idx1, idx2 )

##################################################
# Load CLASH files
##################################################

iFiles = [ os.path.join(data_dir, x) for x in os.listdir(data_dir) if x.endswith('_clash.tsv.gz') ]

for iFile in iFiles:
    clip_data += get_clash_list( iFile )

##################################################
# CRAC data
##################################################

clip_data += get_crac_list( crac_file )

##################################################
# Redundancy Reduce CLIP Data
##################################################

clip_data = redundancy_reduce(clip_data)

##################################################
# Benchmark
###################################################

logger.debug("STEP 6: Benchmarking starbase, clash and crac...")

# benchmark coexpr data
new_coexpr_scores = stringrnautils.benchmark( [x[0] for x in coexpr_data], [x[1] for x in coexpr_data],
                                              [x[2] for x in coexpr_data], [merge_bin(x[5],3,1) for x in coexpr_data],
                                              gold_standard_file_path, discrete=True, fit_name='starbase_coexpr')

coexpr_data = [ x[:5] + [ new_coexpr_scores[idx] ] + x[6:] for idx,x in enumerate(coexpr_data) ]

# benchmark clip data
new_clip_scores = stringrnautils.benchmark( [x[0] for x in clip_data], [x[1] for x in clip_data],
                                            [x[2] for x in clip_data], [merge_bin(x[5],3,1) for x in clip_data],
                                            gold_standard_file_path, discrete=True, fit_name='starbase_clip')


coexpr_data += [ x[:5] + [ new_clip_scores[idx] ] + x[6:] for idx,x in enumerate(clip_data) ]

coexpr_data = redundancy_reduce(coexpr_data)
coexpr_data = [ "\t".join(x[:5] + [ "%0.5f" % x[5]] + x[6:]) + "\n" for x in coexpr_data ]

with open( os.path.join(master_dir, "starbase.tsv"), 'w') as handle:
    handle.write( "".join(coexpr_data)  )

logger.info('Done.' + os.linesep + os.linesep)



