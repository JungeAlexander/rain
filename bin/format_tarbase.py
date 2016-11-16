import gzip, os, argparse, stringrnautils

parser = argparse.ArgumentParser()
parser.add_argument('-data_path',default='data')
parser.add_argument('-rawscore_path', default='rawscore_files')
parser.add_argument('-master_path', default='master_files')
parser.add_argument('-gold_std', default='data/extended_gold_standard.tsv')
args = parser.parse_args()

string_mapper = stringrnautils.get_alias_to_string_mapper(['9606', '10090','7955', '10116', '7227', '6239','3702'], '', '', 10, 'all')
mir_mapper = stringrnautils.get_unique_mir_mapper()

# Retrieve tarbase data
#----------------------
tarbase_file = os.path.join(args.data_path, "Tarbase.6.7.FINAL.mirbase21.download.tsv.gz")
tarbase_file = os.path.join("data", "Tarbase.6.7.FINAL.mirbase21.download.tsv.gz")

if not os.path.exists(tarbase_file):
    os.system( "wget -nv http://rth.dk/~ajunge/Tarbase.6.7.FINAL.mirbase21.download.tsv.gz -O %s" % tarbase_file )

orgn2keep = dict([("Arabidopsis thaliana",'3702'),
                  ("Caenorhabditis elegans",'6239'),
                  ("Danio rerio",'7955'),
                  ("Drosophila melanogaster",'7227'),
                  ("Homo sapiens",'9606'),
                  ("Mus musculus",'10090'),
                  ("Rattus norvegicus",'10116')])

p_idx,r_idx,t_idx,exp_idx,posneg_idx,evidence_idx,direction_idx = 0,2,3,7,8,9,10

def open_file(x):
    return gzip.open(x) if x.endswith(".gz") else open(x)

gold_std_pairs = set([ tuple(x.split("\t")[:3]) for x in open_file(args.gold_std)])

# Define inclusion criterion
def inclusion_criterion(x):
    return  x[exp_idx] == "Luciferase Reporter Assay" and \
            x[t_idx] in orgn2keep and x[evidence_idx] == "DIRECT" and \
            ( (x[direction_idx] == "DOWN" and record[posneg_idx] == "POSITIVE") or record[posneg_idx] == "NEGATIVE")

# Manual correction of miR names
manual_correction = {"hsa-miR-1-3p": "MIMAT0000416",
                     "hsa-miR-203a-3p" :"MIMAT0000264",
                     'hsa-miR-499a-3p':'MIMAT0004772'}

# Define Predictors to validate
predictors = ["miRanda.tsv","miRDB.tsv","starmirdb.tsv", "PITA.tsv","targetscan.tsv","predictions.tsv"]

predictors = [ "%s.gz"%x for x in predictors]

npreds = len(predictors)

# Assemble Tarbase data
failed_mapping = []
data = {}
for record_no, record in enumerate(gzip.open(tarbase_file)):
    if record_no == 0: continue # skip header
    record = record.rstrip().split("\t")
    if not inclusion_criterion(record):
        continue
    taxon = orgn2keep[record[t_idx]]
    failed = [0, 0]
    try:
        if record[r_idx] in manual_correction:
            r_id = mir_mapper[manual_correction[record[r_idx]]]
        else:
            r_id = mir_mapper[record[r_idx]]
    except:
        failed[0]=1
    try:
        p_id = string_mapper[taxon][record[p_idx]]
    except:
        failed[1]=1
    if sum(failed) > 0:
        failed_mapping.append( [record[r_idx], record[p_idx],str(failed[0]),str(failed[1])] )
    else:
        key = (taxon,r_id,p_id)
        if key in gold_std_pairs:
            continue
        if not (record[posneg_idx] == "NEGATIVE" and key in data ): # ignore negative entries if positive entries exist
            data[key] = [ "NA","NA","NA","NA","NA","NA","NA","NA",record[posneg_idx] ]

#print "\n".join(["\t".join(x) for x in failed_mapping])

# Generate data frame with prediction scores

for ifile_no,ifile in enumerate(predictors):
    if ifile_no+1 < npreds:
        ifile = os.path.join(args.rawscore_path,ifile)
        #ifile = os.path.join("rawscore_files",ifile)
    else:
        ifile = os.path.join(args.master_path,ifile)
        #ifile = os.path.join("master_files",ifile)
    for record in open_file(ifile):
        record = record.rstrip().split("\t")
        key = (record[0],record[1],record[2])
        if key in data:
            data[key][ifile_no] = record[5]

output = [ ["Tax","miRNA","Gene"] + predictors + ["Tarbase"] ]
output += [ list(x)+y for x,y in data.iteritems() ]

with open( os.path.join(args.rawscore_path,"Validation_data.tsv"),'w' ) as handle:
    handle.write("".join([ "\t".join(x)+"\n" for x in output ]))
