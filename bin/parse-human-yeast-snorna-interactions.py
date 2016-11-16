import sys,gzip,csv,collections
import stringrnautils
import tarfile
import os.path
import urllib,urllib2
import argparse
import logging
import lxml.html

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

parser = argparse.ArgumentParser(description="""Parses human snoRNAome database and yeast snoRNAs.""")
parser.add_argument('-master_path', help='path to write master files', default = 'master_files')
parser.add_argument('-data_path', help='path to write data files', default = 'data')
parser.add_argument('-id_path', help='path to write data files', default = 'id_dictionaries')
parser.add_argument('--add_rrnas', help='add missing rrna aliases', action='store_true')
parser.add_argument('--add_snornas', help='add missing snorna aliases', action='store_true')


args = parser.parse_args()

id_dict_dir = args.id_path
data_dir = args.data_path
master_dir = args.master_path

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
human_snoRNAome_PATH=data_dir+"/An_updated_human_snoRNAome_Supplementary_Dataset_S1_addnonCanonicalFunc.csv" #interaction file

####### WARNING ##### NOW Parser is fetching the mapping data from public google drive copy but I will try to change that ASAP ######################
human_snoRNAome_URL="http://www.bioinf.uni-leipzig.de/Publications/SUPPLEMENTS/15-065//Supplementary_Dataset_S1_addnonCanonicalFunc.csv"

########### DOwNLOAD data and ID files if not in the data folder
if not os.path.isfile(human_snoRNAome_PATH):
    os.system("wget --no-check-certificate -q -O %s '%s'" % (human_snoRNAome_PATH, human_snoRNAome_URL) )

def get_alias_mapper(ncrna_mapper):
	alias_file_url = "http://rth.dk/resources/rain/data/v1.rna.aliases.complete.txt.gz"
	alias_file = os.path.join(id_dict_dir, "temp.v1.rna.aliases.complete.txt.gz")

	if not os.path.exists(alias_file):
		os.system("wget --no-check-certificate -q -O %s '%s'" % (alias_file, alias_file_url))

	handle = gzip.open(alias_file) if alias_file.endswith(".gz") else open(alias_file)
	tmp_list = [x.strip("\n").split("\t") for x in handle]
	handle.close()
	for l in tmp_list:
		if len(l)>3:
			tax, identifier, alias, source = l
			ncrna_mapper[tax][alias] = identifier
	return ncrna_mapper

def get_rrnas():
	ncrna_mapper = collections.defaultdict(dict)
	alias_file = os.path.join(id_dict_dir, "additional_rrna_alias_file.txt")
	with open(alias_file,'w') as out:
		out.write("9606\t5S_rRNA\t5S\tRAIN\n9606\t5S_rRNA\tRNA5S1\tRAIN\n9606\t5S_rRNA\tRNA5S2\tRAIN\n9606\t5S_rRNA\tRNA5S3\tRAIN\n9606\t5S_rRNA\tRNA5S4\tRAIN\n9606\t5S_rRNA\tRNA5S5\tRAIN\n9606\t5S_rRNA\tRNA5S6\tRAIN\n9606\t5S_rRNA\tRNA5S7\tRAIN\n9606\t5S_rRNA\tRNA5S8\tRAIN\n9606\t5S_rRNA\tRNA5S9\tRAIN\n9606\t5S_rRNA\tRNA5S10\tRAIN\n9606\t5S_rRNA\tRNA5S11\tRAIN\n9606\t5S_rRNA\tRNA5S12\tRAIN\n9606\t5S_rRNA\tRNA5S13\tRAIN\n9606\t5S_rRNA\tRNA5S14\tRAIN\n9606\t5S_rRNA\tRNA5S15\tRAIN\n9606\t5S_rRNA\tRNA5S16\tRAIN\n9606\t5S_rRNA\tRNA5S17\tRAIN\n9606\t5-8S_rRNA\t5-8S\tRAIN\n9606\t5_8S_rRNA\tRNA5-8S1\tRAIN\n9606\t5_8S_rRNA\tRNA5-8S2\tRAIN\n9606\t5_8S_rRNA\tRNA5-8S3\tRAIN\n9606\t5_8S_rRNA\tRNA5-8S4\tRAIN\n9606\t5_8S_rRNA\tRNA5-8S5\tRAIN\n9606\t18S_rRNA\t18S_rRNA\tRAIN\n9606\t18S_rRNA\t18S\tRAIN\n9606\t18S_rRNA\tRNA18S1\tRAIN\n9606\t18S_rRNA\tRNA18S2\tRAIN\n9606\t18S_rRNA\tRNA18S3\tRAIN\n9606\t18S_rRNA\tRNA18S4\tRAIN\n9606\t18S_rRNA\tRNA18S5\tRAIN\n9606\t28S_rRNA\t28S_rRNA\tRAIN\n9606\t28S_rRNA\t28S\tRAIN\n9606\t28S_rRNA\tRNA28S1\tRAIN\n9606\t28S_rRNA\tRNA28S2\tRAIN\n9606\t28S_rRNA\tRNA28S3\tRAIN\n9606\t28S_rRNA\tRNA28S4\tRAIN\n9606\t28S_rRNA\tRNA28S5\tRAIN\n4932\tLSU_rRNA\tLSU_rRNA\tRAIN\n4932\tLSU_rRNA\tLSU\tRAIN\n4932\tSSU_rRNA\tSSU_rRNA\tRAIN\n4932\tSSU_rRNA\tSSU\tRAIN")
	handle = gzip.open(alias_file) if alias_file.endswith(".gz") else open(alias_file)
        tmp_list = [x.strip("\n").split("\t") for x in handle]
        handle.close()
        for l in tmp_list:
                if len(l)>3:
                        tax, identifier, alias, source = l
                        ncrna_mapper[tax][alias] = identifier
        return ncrna_mapper

alias_mapper =  collections.defaultdict(dict) if not args.add_rrnas else get_rrnas()
alias_mapper = get_alias_mapper(alias_mapper)


interactions = set()
snornas_to_add_alias = {}
snornas_to_add_alias["9606"]=set()
snornas_to_add_alias["4932"]=set()

read = False
with open(human_snoRNAome_PATH) as inf:
	reader = csv.reader(inf, delimiter=',', quotechar='"')
	for row in reader:
		if row[0]!='NA' and read:
			snorna = row[0].split('(')[0]
			for i in [20,26,32,38]:
				if row[i]!='NA' and row[i]!='-' and row[i].replace(' ','')!='':
					for inter in row[i].split(";"):
						targetrna = '-'.join(inter.split('-')[:-1])
						if args.add_snornas and snorna not in alias_mapper["9606"].keys():
							alias_mapper["9606"][snorna] = snorna
							snornas_to_add_alias["9606"].add(snorna)
						if snorna in alias_mapper["9606"].keys():
							if targetrna in alias_mapper["9606"].keys():
								interactions.add('9606\t'+alias_mapper["9606"][snorna]+'\t'+alias_mapper["9606"][targetrna]+'\t0\tdatabase\t0.900\tAn updated human snoRNAome\thttps://www.ncbi.nlm.nih.gov/pubmed/27174936')
		if row[0]=='Name':
			read = True


yeast_snorna_url = "http://nar.oxfordjournals.org/content/32/14/4281/T4.expansion.html"
yeast_snorna_page = urllib2.urlopen(yeast_snorna_url)
yeast_snorna_data = yeast_snorna_page.read()

html = lxml.html.fromstring("<table"+yeast_snorna_data.split("</table>")[0].split("<table")[-1]+"</table>")
tbl = []
rows = html.cssselect("tr")
for row in rows:
	tbl.append(list())
	for td in row.cssselect("td"):
		tbl[-1].append(td.text_content())

for row in tbl[2:]:
	targetrna = row[0].split('-')[0]
	if targetrna in alias_mapper["4932"].keys():
		snorna = row[1].split()[0]
		if snorna=="*" or snorna=="" or snorna=="Unknown" or row[6] not in ["1","2","3","4","5"]:
			continue
		if args.add_snornas and snorna not in alias_mapper["4932"].keys():
			alias_mapper["4932"][snorna] = snorna
			snornas_to_add_alias["4932"].add(snorna)
		if snorna in alias_mapper["4932"].keys():
			interactions.add('4932\t'+alias_mapper["4932"][snorna]+'\t'+alias_mapper["4932"][targetrna]+'\t0\tdatabase\t0.900\tsnoRNAdb\thttps://www.ncbi.nlm.nih.gov/pubmed/15306656')

for i in sorted([x for x in interactions]):
	print i
  
if len(snornas_to_add_alias["9606"])>0 or len(snornas_to_add_alias["4932"])>0:
	alias_file = os.path.join(id_dict_dir, "additional_snorna_alias_file.txt")
	with open(alias_file,'w') as out:
		for tax in snornas_to_add_alias.keys():
			for snorna in snornas_to_add_alias[tax]:
				out.write("\t".join([tax,snorna,snorna,"RAIN\n"]))
