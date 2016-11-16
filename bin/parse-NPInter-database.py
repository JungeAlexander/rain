import stringrnautils
import tarfile
import os.path
import argparse
import logging

logger = logging.getLogger(os.path.basename(__file__))
logging.basicConfig(level=logging.INFO)

def is_not_valid_path(arg):
    if os.path.exists(arg):
        parser.error("Given path %s exists. Please provide a non-existing output path or delete the existing file %s manually before running this script." % (arg, arg))
    return arg


def is_valid_path(arg):
    if not os.path.exists(arg):
        parser.error("Given master file path %s does not exist." % arg)
    return arg


parser = argparse.ArgumentParser(description="""Parses ncRNA-RNA ncRNA-Protein interactions from NPINTER database. All mRNAs are mapped to ENSPs using the STRING alias file whenever possible, if there is no conversion, the corresponding interaction is ignored. miRNA identifiers are mapped to those miRBase identifiers used in STRING RNA. mRNA identifier mapping is performed using the stringrnautils module. Others mappings are dependent on the dictionary files those are now downloaded from google drive.""")

parser.add_argument('output_file_path', type=is_not_valid_path, help='The output NPINTER master file is written to this path.')
parser.add_argument('gold_standard_file', type=is_valid_path,
                    help='The gold standard file to benchmark against.')
parser.add_argument('-e', '--exclude_starbase', help='Excludes interactions with PubMed IDs coming from StarBase from the output master file.', action='store_true')
parser.add_argument('-i', '--include_starbase', help='Includes interactions with PubMed IDs coming from StarBase in the output master file.', action='store_true')
parser.add_argument('-master_path', help='path to write master files', default = 'master_files')
parser.add_argument('-data_path', help='path to write data files', default = 'data')
parser.add_argument('-id_path', help='path to write data files', default = 'id_dictionaries')

args = parser.parse_args()

output_file_path = args.output_file_path
include_starbase = args.include_starbase
exclude_starbase = args.exclude_starbase
id_dict_dir = args.id_path
data_dir = args.data_path
master_dir = args.master_path
gold_standard_file_path = args.gold_standard_file

if (not include_starbase and not exclude_starbase) or (include_starbase and exclude_starbase):
    raise IOError("Please select exactly one of the following parameters: --exclude_starbase or --include_starbase.")

##########################
# Generate Directory Structure
##########################

if not os.path.exists(id_dict_dir):
    os.mkdir(id_dict_dir)

if not os.path.exists(data_dir):
    os.mkdir(data_dir)

if not os.path.exists(master_dir):
    os.mkdir(master_dir)
##############################33


#FILE PATHS
NP_PATH=data_dir+"/interaction_NPInter[v2.0].txt" #interaction file
BIOMART_DIC_PATH=id_dict_dir+"/NPINTERmapfiles.zip"
NONCODE_DIC_PATH=id_dict_dir+"/NONCODEID_conversions.tsv"
GENENAME_DIC_PATH=id_dict_dir+"/ncRNAaliasfile.tsv.gz"

COMMENT = ''

if not os.path.isfile(NP_PATH):
    tfile = tarfile.open("data/interaction_NPInter[v2.0].txt.tar.gz", 'r:gz')
    tfile.extractall('data/')

# Organism names in NPINTER and their STRING ID conversions
organismIdMap = {
    'Caenorhabditis elegans': '6239',
    'Danio rerio':'7955',
    'Drosophila melanogaster':'7227',
    'Homo sapiens': '9606',
    'Mus musculus':'10090',
    'Oryctolagus cuniculus':'9986',
    'Saccharomyces cerevisiae':'4932',
    'Saccharomyces cerevisiae S288c':'4932'

    # FOLLOWING ORGANISMS DONT HAVE BIOMART Dictionaries
    #'Salmonella typhimurium':'99287'
    #'Candida albicans SC5314':'5476',
    #'Cricetulus griseus':'',
    #'Dictyostelium discoideum AX4':'44689',
    #'Escherichia coli':'',
    #'Kaposi sarcoma-associated herpesvirus (KSHV)':'37296',
    #'Scarabaeus aegyptiorum':'',
    #'Staphylococcus aureus':'',
    #'Arabidopsis thaliana':'3702',
    #'Bacillus subtilis': '224308',
}

organismAbbrevations={'ce':'Caenorhabditis elegans',
                      'dr':'Danio rerio',
                      'dm':'Drosophila melanogaster',
                      'hs':'Homo sapiens',
                      'mm':'Mus musculus',
                      'oc':'Oryctolagus cuniculus',
                      'sc':'Saccharomyces cerevisiae'}

excludedPubMedIds = set()
if exclude_starbase:
    logger.info('STARBASE EXCLUDED')
    excludedPubMedIds =  stringrnautils.starbase_exp_pmids()

#21037263 is the StarBase PubMed ID. 184254 of 201107 interactions in NPINTER are from there.
#These are the other popular interaction sources. 21358643 (1320); 22473208 (2663); 20371350 (2295); 22081015 (1048).


NPINTERorganisms=set()
interactions={}
DescriptionDic = {}

########### Get the number of interactions for all pubmed ids in the database
def getScoresForPubMedIds():
    npinter_f = open(NP_PATH,'r')
    names = {}
    for line in npinter_f:
        interID,ncID,ncType,ncIdentifier,ncName,PartnerID,prType,prIdentifier,prName,interDescription,experiment,pubmed,organism,tag,interClass,interLevel = line.rstrip().split('\t')
        if names.has_key(pubmed):
            names[pubmed]=names[pubmed]+1
        else:
            names[pubmed]=1
    npinter_f.close()
    return names


########## STRING dic for ENSPlike conversions
def getSTRINGdic(specie):					#ENSG ENSP conversion for RefSeq NM_ mRNAs
    STRING_dic = stringrnautils.get_alias_to_string_mapper(organisms=organismIdMap[specie],filter_string_alias='', filter_string_id='')
    return STRING_dic

######### Parse with Specie Name ################################################
UP_dic = stringrnautils.getUniProtDic(BIOMART_DIC_PATH)
NM_dic = stringrnautils.getRefSeqNMdic(BIOMART_DIC_PATH)
NR_dic = stringrnautils.getRefSeqNRdic(BIOMART_DIC_PATH,GENENAME_DIC_PATH)
NONCODE_dic = stringrnautils.getNONCODEdic(NONCODE_DIC_PATH,BIOMART_DIC_PATH,GENENAME_DIC_PATH)
MB_dic = stringrnautils.get_unique_mir_mapper()

def ParseNPINTER(specie = 'Homo sapiens'):
    #for Homo Sapiens if specie is not chosen
    #RNAs or Proteins in NPInter has identifiers from NonCode(NR_, ENST), RefSeq(NM_), MirBase(miR) and UniProt.
    unmapped = {}
    mapped = {}
    unc=0
    mac=0
    UniProtdic={}
    if UP_dic.has_key(organismIdMap[specie]):
        UniProtdic = UP_dic[organismIdMap[specie]]
    RefSeqdic={}
    if NM_dic.has_key(organismIdMap[specie]):
        RefSeqdic = NM_dic[organismIdMap[specie]]
    NC_dic={}
    if NONCODE_dic.has_key(organismIdMap[specie]):
        NC_dic = NONCODE_dic[organismIdMap[specie]]
    logger.info('Parsing NPINTER for ' + specie+ '. UniProt dict size:' + str(len(UniProtdic)) + ' RefSeq dict size:' +
                                                        str(len(RefSeqdic)) + ' NONCODE dict size: '+ str(len(NC_dic)))
    npinter_f = open(NP_PATH,'r')
    npinter_f.readline()
    experiments = set()
    for line in npinter_f:
        interID,ncID,ncType,ncIdentifier,ncName,PartnerID,prType,prIdentifier,prName,interDescription,experiment,pubmed,organism,tag,interClass,interLevel = line.rstrip().split('\t')
        if organism.lower().startswith(specie.lower()) and (interLevel.startswith('RNA-RNA') or interLevel.startswith('RNA-Protein')):
            source='NA'
            target='NA'
            interDescription=interDescription.rstrip()
            if ncType=='miRBase':
                if MB_dic.has_key(ncIdentifier):
                    source = MB_dic[ncIdentifier]
                #else:
                #	sys.stdout.write(ncType+'\t'+ncIdentifier+'\t'+ncName+'\n')
            elif ncType=='RefSeq':
                ID = ncIdentifier.split(".")[0]
                if RefSeqdic.has_key(ID):
                    source = RefSeqdic[ID]
                elif RefSeqdic.has_key(ncIdentifier):
                    source = RefSeqdic[ncIdentifier]
                #else:
                #	sys.stdout.write(ncType+'\t'+ncIdentifier+'\t'+ncName+'\n')
            elif ncType=='UniProt':
                if UniProtdic.has_key(ncName):
                    source = UniProtdic[ncName]
                elif UniProtdic.has_key(ncIdentifier):
                    source = UniProtdic[ncIdentifier]
                    #else:
                    #	sys.stdout.write(ncType+'\t'+ncIdentifier+'\t'+ncName+'\n')
            elif ncType=='NONCODE':
                if NC_dic.has_key(ncIdentifier):
                    source = NC_dic[ncIdentifier]
                    #else:
                    #	sys.stdout.write(ncType+'\t'+ncIdentifier+'\t'+ncName+'\n')
            else:
                logger.debug(ncType)
            if prType=='miRBase':
                if MB_dic.has_key(prIdentifier):
                    target = MB_dic[prIdentifier]
                    #else:
                    #	sys.stdout.write(prType+'\t'+prIdentifier+'\t'+prName+'\n')
            elif prType=='RefSeq':
                ID = prIdentifier.split(".")[0]
                if RefSeqdic.has_key(prIdentifier):
                    target = RefSeqdic[prIdentifier]
                elif RefSeqdic.has_key(ID):
                    target = RefSeqdic[ID]
                #else:
                #	sys.stdout.write(prType+'\t'+prIdentifier+'\t'+prName+'\n')
                elif prType=='UniProt':
                    if UniProtdic.has_key(prName):
                        target = UniProtdic[prName]
                    elif UniProtdic.has_key(prIdentifier):
                        target = UniProtdic[prIdentifier]
                        #else:
                        #	sys.stdout.write(prType+'\t'+prIdentifier+'\t'+prName+'\n')
                elif prType=='NONCODE':
                    if NC_dic.has_key(prIdentifier):
                        target = NC_dic[prIdentifier]
                        #else:
                        #	sys.stdout.write(prType+'\t'+prIdentifier+'\t'+prName+'\n')
                else:
                    logger.debug(prType)
                if source!='NA' and target!='NA':
                    if mapped.has_key(ncType+' '+prType):
                        mapped[ncType+' '+prType] += 1
                    else:
                        mapped[ncType+' '+prType] = 1
                    mac+=1
                    if interactions.has_key(source+target) or interactions.has_key(target+source):
                        if interactions.has_key(source+target):
                            (orgid,dicid1,dicid2,pubmeds,score,descriptions)=interactions[source+target]
                        else:
                            (orgid,dicid1,dicid2,pubmeds,score,descriptions)=interactions[target+source]
                        pubmeds.add(pubmed)
                        score = score + (1.0/PubMedCounts[pubmed])
                        if DescriptionDic.has_key(interDescription):
                            descriptions.add(DescriptionDic[interDescription])
                        else:
                            descriptions.add('Unknown')
                            #sys.stdout.write(interDescription+'\n\n')
                            interactions[dicid1+dicid2]=(orgid,dicid1,dicid2,pubmeds,score,descriptions)
                    else:
                        score = 1.0/PubMedCounts[pubmed]
                        if DescriptionDic.has_key(interDescription):
                            interactions[source+target]=(organismIdMap[specie],source,target,
                                                         {pubmed},score, {DescriptionDic[interDescription]})
                        else:
                            interactions[source+target]=(organismIdMap[specie],source,target,
                                                         {pubmed},score, {'Unknown'})
                            #sys.stdout.write(interDescription+'\n\n')
        else:
            if unmapped.has_key(ncType+' '+prType):
                unmapped[ncType+' '+prType] += 1
            else:
                unmapped[ncType+' '+prType] = 1
                unc+=1
    return (unmapped,mapped,unc,mac)

################# Benchmark
def BenchmarkWithMovingAvg():	#interaction list obtained by parseNPINTER()
    rnaIds=[]
    prIds=[]
    scores=[]
    for (orgID,id1,id2,pubmeds,rawscore,descs) in interactions.values():
        if id1.startswith('hsa') and id2.startswith('ENSP'):
            rnaIds.append(id1)
            prIds.append(id2)
            scores.append(float(rawscore))
        elif id2.startswith('hsa') and id1.startswith('ENSP'):
            rnaIds.append(id2)
            prIds.append(id1)
            scores.append(float(rawscore))

    stringrnautils.benchmark(['9606'] * len(rnaIds), rnaIds, prIds, scores, gold_standard_file_path, window_size=10)

############# Benchmarking is done by two bins Single Source or Multiple Source
def BenchmarkWithBins():
    rnaIds=[]
    prIds=[]
    bins=[]
    organisms=[]
    for (orgID,id1,id2,pubmeds,rawscore,descs) in interactions.values():
        rnaIds.append(id1)
        prIds.append(id2)
        #bins.append(str(len(pubmeds))+' pubmed ref')
        if len(pubmeds)==1:
            bins.append('Single Source')
        else:
            bins.append('Multiple Source')
        organisms.append(orgID)

    scoreslist = stringrnautils.discrete_benchmark(organisms, rnaIds, prIds, bins, gold_standard_file_path)
    scores = {}
    n=0
    while n < len(rnaIds):
        score = 0
        while (n+1)<len(rnaIds) and rnaIds[n]==rnaIds[n+1] and prIds[n]==prIds[n+1]:
            score += scoreslist[n]
            n+=1
        score += scoreslist[n]
        scores[rnaIds[n]+prIds[n]]=score
        n+=1

    return scores

########### Parsing

logger.info("Parsing NPInter.")

#DescriptionDic = getDescriptionDic()
PubMedCounts = getScoresForPubMedIds()

for ab in organismAbbrevations.keys():
    specie = organismAbbrevations[ab]
    NPINTERorganisms.add(organismIdMap[specie])
    un,mapped,unc,mac = ParseNPINTER(specie)
    logger.info(str(specie) + ': Unmapped Interactions'+str(un))
    logger.info(str(specie) + ': Mapped Interactions'+str(mapped))

#Get Scores
scores = BenchmarkWithBins()

#Write interactions excluding Starbase and the other banned pubmed ids
outf=open(output_file_path,'w')
count=0
for (orgID,id1,id2,pubmeds,rawscore,descs) in interactions.values():
    url="http://www.ncbi.nlm.nih.gov/pubmed/"
    valid = False
    for pubmed in pubmeds:
        if int(pubmed) not in excludedPubMedIds:
            url+=pubmed+','
            valid=True
    if valid:
        count+=1
        outf.write(orgID+'\t'+id1+'\t'+id2+'\t'+'0'+'\t'+'Experiments'+'\t'+str(scores[id1+id2])+'\t'+'NPInter'+'\t'+url.rstrip(',')+ '\t' + COMMENT +'\n')
outf.close()

logger.info('Adding ' + str(count) + ' interactions in total')
logger.info('Done.' + os.linesep + os.linesep)
