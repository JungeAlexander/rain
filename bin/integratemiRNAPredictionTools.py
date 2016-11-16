import gzip
import os
import collections
import argparse
import stringrnautils
import logging
logger = logging.getLogger(os.path.basename(__file__))
logging.basicConfig(level=logging.INFO)


def is_valid_path(arg):
    if not os.path.exists(arg):
        parser.error("Given file/dir path %s does not exist." % arg)
    return arg

parser = argparse.ArgumentParser(description="""
                                    Parses and integrates the miRNA-target prediction tools TargetScan, miRanda,
                                    PicTar, and STarMirDB.
                                    Interactions predicted by these tools are collected in a single master file
                                    with each interaction being scored using functions in bin/stringrnautils.py.
                                    """)

parser.add_argument('gold_standard_file', type=is_valid_path,
                    help='The gold standard file to benchmark against.')
parser.add_argument('-master_path',
                    help='Output master file is written to the given folder.',
                    default='master_files')
parser.add_argument('-rawscore_path',
                    help='Output raw scores in master file format is written to the given folder.',
                    default='rawscore_files')
parser.add_argument('-data_path',
                    help='Raw data files are written to the given folder.',
                    default='data/')
parser.add_argument('-run_all',
                    help='Also read and write raw scores of predictors that are excluded from RAIN',
                    action='store_true')
parser.add_argument('-species_norm',
                    help='Normalize across species prior to benchmarking',
                    action='store_true')
args = parser.parse_args()

LOCAL_DATA_PATH = args.data_path
MASTER_FILE_DIR = args.master_path
RAWSCORE_FILE_DIR = args.rawscore_path
gold_standard_file_path = args.gold_standard_file

if not os.path.exists(RAWSCORE_FILE_DIR):
    os.mkdir(RAWSCORE_FILE_DIR)

PREDICTION_MASTER_FILE_NAME = 'prediction' # evidence name in master files


#################
# Methods
#################

def species_covered(maser_file):
    return {int(tab.split('\t', 1)[0]) for tab in open(maser_file)}


def write2master( filename, method, tax_ids,rna_ids,protein_ids,scores,master_dir=MASTER_FILE_DIR ):
    with open( os.path.join(master_dir, filename), 'w' ) as fw:
        for organism_id, rna_id, protein_id, new_score in zip(tax_ids, rna_ids, protein_ids, scores):
            fw.write('\t'.join((organism_id, rna_id, protein_id, '1',
                                PREDICTION_MASTER_FILE_NAME,
                                str(new_score), method, '', ''))+'\n')

def starmirdb_subfunc(filenames, mir_mapper, string_mapper, ofile_base):
    """
    function: read_starmirdb - Read and process StaRMIRDB predictions

    Arguments       Type    Description
    --------------------------------
    filenames       list    List of files containing the starmirdb data
    mir_mapper      dict    Dictionary for mapping of ids to STRING-RNA miR ids, wanted_id = hash[id]
    string_mapper   dict    Dictionary for mapping of protein and gene ids to STRING ids, wanted_id = hash[tax][id]
    ofile_base      string  Base name of the out files

    Author: Christian Garde (garde@cbs.dtu.dk), Workman Lab, DTU Systems Biology
    Last revision April 29, 2015

    """

    #------------
    # Read Raw Data
    #------------
    idx1, idx2, idx_score = (0, 1, 41)
    data = collections.defaultdict(dict)
    verbose = False # major speed cost, will collect mapping statistics
    miRs_mapped,miRs_not_mapped, proteins_mapped, proteins_not_mapped = [],[],[],[]

    for filename in filenames:
        filename = os.path.join(LOCAL_DATA_PATH, filename)
        tax_id = filename.split("/")[-1].split("_")[0]
        handle = gzip.open(filename) if filename.endswith(".gz") else open(filename)
        for idx, record in enumerate(handle):
            if idx == 0: # header
                continue
            record = record.split("\t")
            mir_id = 'mmu-' + record[idx2] if tax_id == '10090' else record[idx2]
            if not mir_id in mir_mapper:
                if verbose:
                    miRs_not_mapped.append(mir_id)
                continue
            if verbose:
                miRs_mapped.append(mir_id)
            id1 = mir_mapper[mir_id]
            protein_id = record[idx1]
            if not protein_id in string_mapper[tax_id]:
                if verbose:
                    proteins_not_mapped.append(tax_id + protein_id)
                continue
            if verbose:
                proteins_mapped.append(tax_id + protein_id)
            id2 = string_mapper[tax_id][protein_id]
            key = (id1, id2)
            score = float(record[idx_score])
            if key in data[tax_id]:
                data[tax_id][key].append( score )
            else:
                data[tax_id][key] = [ score ]

    if verbose:
        miRs_mapped = set(miRs_mapped)
        miRs_not_mapped = set(miRs_not_mapped)
        proteins_mapped = set(proteins_mapped)
        proteins_not_mapped = set( proteins_not_mapped )
        #
        logging.info('STarMirDB mapping statistics:')
        logging.info('-----------------------------')
        logging.info('miRs mapped to STRING-RNA IDs using miR mapping in bin/stringrnautils.py:%i' %
                         len(miRs_mapped))
        logging.info('miRs NOT mapped to STRING-RNA IDs using miR mapping in bin/stringrnautils.py:%i' %
                         len(miRs_not_mapped))
        miRs_not_mapped_indented = ['\t' + m for m in miRs_not_mapped]
        logging.info('miRs not mapped were:\n%s' % '\n'.join(miRs_not_mapped_indented))
        logging.info('')

        logging.info('Proteins mapped to STRING-RNA IDs using protein mapping in bin/stringrnautils.py:%i' %
                         len(proteins_mapped))
        logging.info('Proteins NOT mapped to STRING-RNA IDs using protein mapping in bin/stringrnautils.py:%i' %
                         len(proteins_not_mapped))
        proteins_not_mapped_indented = ['\t' + m for m in proteins_not_mapped]
        logging.info('Proteins not mapped were:\n%s' % '\n'.join(proteins_not_mapped_indented))
        logging.info('')

    #------------
    # Remove duplicates
    #------------
    for tax in data.keys():
        stringrnautils.reduce_dict_scores( data[tax], "max" )

    #------------
    # normalize across species
    #------------
    if args.species_norm:
        data = stringrnautils.qq_correct( data, "%s_normalization.pdf"%ofile_base)

    tax_ids,rna_ids,protein_ids,scores = zip(*[ (tax, mol_ids[0], mol_ids[1], score)
                                         for tax,mol_dict in data.iteritems()
                                         for mol_ids,score in mol_dict.iteritems() ])
    #------------
    # Benchmark data
    #------------
    new_scores = stringrnautils.benchmark(tax_ids, rna_ids, protein_ids, scores, gold_standard_file_path,
                                          window_size=200, fit_name=ofile_base)
    #------------
    # Write data to master file
    #------------
    write2master( '%s.tsv'%ofile_base, 'STarMirDB', tax_ids, rna_ids, protein_ids, new_scores, MASTER_FILE_DIR )
    write2master( '%s.tsv'%ofile_base, 'STarMirDB', tax_ids, rna_ids, protein_ids, scores, RAWSCORE_FILE_DIR )

def read_starmirdb(mir_mapper, string_mapper):
    """
    function: read_starmirdb - Read and process StaRMIRDB predictions

    Arguments       Type    Description
    --------------------------------
    mir_mapper      dict    Dictionary for mapping of ids to STRING-RNA miR ids, wanted_id = hash[id]
    string_mapper   dict    Dictionary for mapping of protein and gene ids to STRING ids, wanted_id = hash[tax][id]

    Author: Christian Garde (garde@cbs.dtu.dk), Workman Lab, DTU Systems Biology
    Last revision April 29, 2015

    """

    #------------
    # Retrieve data from RTH page
    #------------
    starmirdb_file = os.path.join(LOCAL_DATA_PATH, 'starmirdb.tar.gz')
    os.system("tar -zxf %s -C %s"%(starmirdb_file, LOCAL_DATA_PATH))

    CDS_files = ["9606_starmirdb_CDS.tsv.gz","10090_starmirdb_CDS.tsv.gz"]
    UT3_files = ["9606_starmirdb_UTR3.tsv.gz","10090_starmirdb_UTR3.tsv.gz"]
    UTR5_files = ["9606_starmirdb_UTR5.tsv.gz","10090_starmirdb_UTR5.tsv.gz"]

    #starmirdb_subfunc(CDS_files, mir_mapper, string_mapper, "starmirdb_CDS")
    starmirdb_subfunc(UT3_files, mir_mapper, string_mapper, "starmirdb")
    #starmirdb_subfunc(UTR5_files, mir_mapper, string_mapper, "starmirdb_UTR5")

def read_predictions( pred_file, gene2ensembl, mir_mapper, string_mapper, tax_idx=0, mir_idx=1,target_idx=2,score_idx=3,
                      increasing=True, window_size=50, name="miRDB", ignore_fraction=0.0, has_header=True, do_benchmark=True ):
    """
    function: read_pictar - read and process PicTar predictions

    Arguments       Type    Description
    --------------------------------
    ifile           str     Name of file to retrieve
    gene2ensembl    dict    Dictionary for mapping gene names to ensemble id ensembl_id = hash[gene_name], if len(gene2ensembl)=0 it's ignored
    mir_mapper      dict    Dictionary for mapping of ids to STRING-RNA miR ids, wanted_id = hash[id]
    string_mapper   dict    Dictionary for mapping of protein and gene ids to STRING ids, wanted_id = hash[tax][id]
    tax_idx         int     Idx in tab separated file of the taxonomy identifier
    mir_id          int     Idx in tab separated file of the miRNA identifier
    target_idx      int     Idx in tab separated file of the mRNA target identifier
    score_idx       int     Idx in tab separated file of the prediction score
    increasing      bool    Whether scoring scheme defines high as better or low as better
    window_size     int     Size of sliding window during benchmarking
    name            str     Name used for the masterfile and benchmarking plots
    ignore_fraction float   Burn in before fitting
    has_header      bool    Whether the first line in the file should be interpreted as a header
    do_benchmark    bool    Whether to benchmark or simply write out raw scores

    Last revision November 27, 2015
    Christian Garde

    """
    # Retrieve data
    #---------
    data_file = os.path.join(LOCAL_DATA_PATH,pred_file)

    def open_file(x):
        return gzip.open(x) if x.endswith(".gz") else open(x)
    #
    # Load data
    #---------
    data = collections.defaultdict(dict)
    for idx, record in enumerate(open_file(data_file)):
        if idx==0 and has_header: continue # skip header
        record = record.rstrip('\n').split('\t')
        tax,miR,gene,score = record[tax_idx],record[mir_idx],record[target_idx],record[score_idx]
        try:
            if False:#len(gene2ensembl)>0:
                gene = gene2ensembl[tax][gene]
            key = (mir_mapper[miR], string_mapper[tax][gene])
            score = float(score)
        except:
            continue
        if key in data[tax]:
            data[tax][key].append( score )
        else:
            data[tax][key] = [ score ]
    #
    # Redundancy Reduce
    #---------
    for tax in data.keys():
        stringrnautils.reduce_dict_scores( data[tax], "max" if increasing else "min" )

    # Species Normalization
    #---------
    #if args.species_norm:
    #    data = stringrnautils.qq_correct( data, '%s_normalization.pdf' % name)
    #
    # Benchmark
    #--------
    tax_ids,rna_ids,protein_ids,scores = zip(*[ (tax, mol_ids[0], mol_ids[1], score)
                                         for tax,mol_dict in data.iteritems()
                                         for mol_ids,score in mol_dict.iteritems() ])
    # Write raw scores to file
    write2master( '%s.tsv' % name, name, tax_ids, rna_ids, protein_ids, scores, RAWSCORE_FILE_DIR )
    #
    if do_benchmark:
        new_scores  = stringrnautils.benchmark(tax_ids, rna_ids, protein_ids,
                                               scores if increasing else [-1*x for x in scores ],
                                               gold_standard_file_path,increases=True,
                                               window_size=window_size, fit_name=name,
                                               ignore_fraction=ignore_fraction)
        # Write to masterfile
        write2master( '%s.tsv' % name, name, tax_ids, rna_ids, protein_ids, new_scores, MASTER_FILE_DIR )


def integrate_all_prediction_tools():
    # Define dictionaries
    #--------------------
    gene2ensembl = stringrnautils.map_gene_2_enemble(os.path.join(LOCAL_DATA_PATH, 'gene2ensembl.gz'))
    stringrnautils.integrate_NM_dictionary(gene2ensembl)

    mir_mapper = stringrnautils.get_unique_mir_mapper()
    string_mapper = stringrnautils.get_alias_to_string_mapper(['9606', '10090','7955', '10116', '7227', '6239','3702'], '', '', 10, 'all')

    # Read data and benchmark
    #--------------------------

    # starmirdb - may decide to exclude this one
    read_starmirdb( mir_mapper, string_mapper)

    # miRanda
    read_predictions( "miRanda_v3.3a.tsv.gz", {}, mir_mapper, string_mapper,
                      tax_idx=0, mir_idx=1,target_idx=2,score_idx=4,
                      increasing=False, window_size=1000, name="miRanda",
                      ignore_fraction=0.7, has_header=True )

    # miRDB
    read_predictions( "miRDB_v5.0.tsv.gz", gene2ensembl, mir_mapper, string_mapper,
                      tax_idx=0, mir_idx=1,target_idx=2,score_idx=3,
                      increasing=True, window_size=75, name="miRDB",
                      ignore_fraction=0.0, has_header=True )

    # PITA
    read_predictions( "PITA.tsv.gz", {}, mir_mapper, string_mapper,
                      tax_idx=0, mir_idx=2,target_idx=1,score_idx=4,
                      increasing=False, window_size=500, name="PITA",
                      ignore_fraction=0.0, has_header=True )

    # RNA22 - excluded due to poor performance
    if args.run_all:
        read_predictions( "RNA22.tsv.gz", {}, mir_mapper, string_mapper,
                          tax_idx=0, mir_idx=1,target_idx=2,score_idx=3,
                          increasing=True, window_size=50, name="RNA22",
                          ignore_fraction=0.2, has_header=True,do_benchmark=False)

    # RNAhybrid - excluded due to poor performance
    if args.run_all:
        read_predictions( "RNAhybrid_seed.tsv.gz", {}, mir_mapper, string_mapper,
                          tax_idx=0, mir_idx=1,target_idx=2,score_idx=3,
                          increasing=False, window_size=50, name="RNAhybrid_seed",
                          ignore_fraction=0.2, has_header=False, do_benchmark=False)

    # Targetscan
    read_predictions( "targetscan.mammals.tsv.gz", {}, mir_mapper, string_mapper,
                      tax_idx=0, mir_idx=1,target_idx=2,score_idx=4,
                      increasing=False, window_size=50, name="targetscan",
                      ignore_fraction=0.50, has_header=True )

    # integrate prediction tools
    #--------------------
    prediction_tools = ('starmirdb', 'miRanda', 'targetscan', 'miRDB', 'PITA')
    organism_to_tool = {}
    for tool in prediction_tools:
        organisms = species_covered(os.path.join(MASTER_FILE_DIR,'{0}.tsv'.format(tool)))
        for organism in organisms:
            organism_to_tool.setdefault(organism, []).append(tool)

    tool_combinations = set()
    tool_combinations_to_species = {}
    for organism, tools in list(organism_to_tool.items()):
        tools = '_and_'.join(sorted(tools))
        tool_combinations.add(tools)
        organism_to_tool[organism] = tools
        tool_combinations_to_species.setdefault(tools, set()).add(organism)

    tool_parameters = {
    "PITA_and_miRDB_and_miRanda_and_targetscan":{
        'negative_evidence' : False,
        'rebenchmark_everything' : True,
        'ignore_fraction' : 0.0,
        'window_size' : 110,
        'unlink_master_files' : False
    },
    "PITA_and_miRanda": {
        'negative_evidence' : False,
        'rebenchmark_everything' : True,
        'ignore_fraction' : 0.0,
        'window_size' : 200,
        'unlink_master_files' : False
    }
    }

    default_tool_parameters = {
        'negative_evidence' : False,
        'rebenchmark_everything' : True,
        'ignore_fraction' : 0.60,
        'window_size' : 75,
        'unlink_master_files' : False
    }

    # generate organism specific callibration curves
    predictions_master_file = open(os.path.join(MASTER_FILE_DIR, 'predictions.tsv'), 'w')
    new_master_files = ['{0}.tsv'.format(p) for p in prediction_tools]

    for tool_combination in tool_combinations:
        source_master_files = ('{0}.tsv'.format(t) for t in tool_combination.split('_and_'))
        destination_name = 'predictions_subset_{0}'.format(tool_combination)
        destination_master_file = 'predictions_subset_{0}.tsv'.format(tool_combination)

        parameters = default_tool_parameters.copy()
        if tool_combination in tool_parameters:
            parameters.update(tool_parameters[tool_combination])

        new_master_files.append(destination_master_file)
        stringrnautils.combine_masterfiles(source_master_files, destination_master_file,
                                           gold_standard_file_path, destination_name,
                                           **parameters)

        # generate/append relevant species to predictions.tsv
        species = tool_combinations_to_species[tool_combination]
        for line in open(os.path.join(MASTER_FILE_DIR, destination_master_file)):
            if int(line.split('\t', 1)[0]) in species:
                predictions_master_file.write(line)

    # delete all the tmp master files
    for master_file in new_master_files:
        os.unlink(os.path.join(MASTER_FILE_DIR, master_file))


    # stringrnautils.combine_masterfiles(('starmirdb.tsv',
    #                                     'miRanda.tsv',
    #                                     'targetscan.tsv',
    #                                     'miRDB.tsv',
    #                                     'PITA.tsv'),
    #                                    'predictions.tsv', gold_standard_file_path, 'predictions',75,
    #                                    negative_evidence=False,rebenchmark_everything=True,
    #                                    ignore_fraction=0.60)

if __name__ == '__main__':

    ##################
    # EXECUTE CODE
    ##################
    logger.info('Parsing prediction tools')
    integrate_all_prediction_tools()
    logger.info('Done.' + os.linesep + os.linesep)
