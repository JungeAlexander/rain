import argparse
import gzip
import hashlib
import logging
import os
import requests
import sys
import urllib

logger = logging.getLogger(os.path.basename(__file__))
logging.basicConfig(level=logging.INFO)
logging.getLogger('requests').setLevel(logging.WARNING)  # to suppress HTTP connection info

__author__ = 'Alexander Junge (alexander.junge@gmail.com), RTH, University of Copenhagen'

parser = argparse.ArgumentParser(description='''
Checks if all input data files needed by the RAIN pipeline exist locally
and downloads missing files if necessary. If downloading any of the non-existing
files fails, the present script exists with a non-zero exit status.
 ''')
parser.add_argument('-f', '--force', action='store_true',
                    help='If set, all data files are downloaded, also files existing locally.')
parser.add_argument('-s', '--skip', action='store_true',
                    help='DANGER ZONE: If set, files to be downloaded that do not exist are skipped. Careful with '
                         'setting this as it may lead to an altered/incomplete RAIN pipeline output.')
args = parser.parse_args()

force_dl = args.force
skip_unavailable_dl = args.skip


def check_and_download(path, url, force_download=False, skip_unavailable=False,
                       md5sum=None):
    """
    Check if a file exists and download if necessary.

    :param path: the location of the file to be checked
    :param url: the URL the file is to be downloaded from
    :param force_download: Default: False. Set to True to force file download irrespective of file existence.
    :param skip_unavailable: Default: False. Set to True to skip files that are not available for download or
                whose local md5 sum does not match
    :param md5sum: string; default: None. The expected md5 hash sum of the file to download
                If the md5 of the file at path does not match this expected sum,
                the script fails.
    """
    if not os.path.exists(path) or force_download:
        logger.info('Downloading dependency {} from URL: {}'.format(path, url))
        if url.startswith('http'):
            r = requests.get(url, stream=True)
            if r.status_code == 200:
                with open(path, 'wb') as f:
                    for chunk in r.iter_content(1024):
                        f.write(chunk)
            else:
                if not skip_unavailable:
                    logger.error('Exiting! Could not download dependency {} from URL: {}'.format(path, url))
                    sys.exit(1)
                else:
                    logger.warning('Could not download dependency {} from URL: {}'.format(path, url))
        else:
            try:
                urllib.urlretrieve(url, path)
            except IOError as e:
                if not skip_unavailable:
                    logger.error('Exiting! Could not download dependency {} from URL: {}'.format(path, url))
                    raise e
                else:
                    logger.warning('Could not download dependency {} from URL: {}'.format(path, url))
    else:
        logger.info('Not downloading dependency {} as it was found locally.'.format(path))

    if md5sum is not None:
        computed_md5 = hashlib.md5(open(path, 'rb').read()).hexdigest()
        if computed_md5 != md5sum and not skip_unavailable:
            logger.error('Exiting! MD5 hash {} for local file {} did not match '
                         'expected hash {}'.format(computed_md5, path, md5sum))
            logger.error('Please delete local file {} an restart this '
                         'script to download the file again.'.format(path))
            sys.exit(1)

id_dict_dir = 'id_dictionaries'
data_dir = 'data'

miRTarBase_version = '6.1'
if miRTarBase_version == '4.5':
    miRTarBase_URL = 'http://mirtarbase.mbc.nctu.edu.tw/cache/download/' + miRTarBase_version + '/miRTarBase_MTI.xls'
    miRTarBase_file = os.path.join(data_dir, 'miRTarBase_' + miRTarBase_version + '_MTI.xls')
elif miRTarBase_version == '6.1':
    miRTarBase_URL = 'http://mirtarbase.mbc.nctu.edu.tw/cache/download/6.1/miRTarBase_MTI.xlsx'
    miRTarBase_file = os.path.join(data_dir, 'miRTarBase_' + miRTarBase_version + '_MTI.xlsx')
else:
    raise ValueError('miRTarBase version %s is not supported.' % miRTarBase_version)

NPInter_URL = 'http://www.bioinfo.org/NPInter/datadownload/interaction_NPInter%5Bv2.0%5D.txt.tar.gz'
NPInter_file = os.path.join(data_dir, 'interaction_NPInter[v2.0].txt.tar.gz')

mirbase_version = 20
mirbase_base_url = 'ftp://mirbase.org/pub/mirbase/%i'
base_url = mirbase_base_url % mirbase_version
mirbase_aliases_txt_base_path = os.path.join(id_dict_dir, 'miRBase_%i_aliases.txt.gz')
mirbase_aliases_txt_path = mirbase_aliases_txt_base_path % mirbase_version
mirbase_aliases_txt_url = '%s/aliases.txt.gz' % base_url
mirbase_mirna_dat_base_path = os.path.join(id_dict_dir, 'miRBase_%i_miRNA.dat.gz')
mirbase_mirna_dat_path = mirbase_mirna_dat_base_path % mirbase_version
mirbase_mirna_dat_url = '%s/miRNA.dat.gz' % base_url
mir_mature_file_path = os.path.join(id_dict_dir, 'mature.fa.gz')
mir_mature_file_url = '%s/mature.fa.gz' % base_url
mir_fam_url = '%s/miFam.dat.gz' % base_url
mir_fam_path = os.path.join(data_dir, 'miRBase_%i_miFam.dat.gz' % mirbase_version)

text_mining_urls_tax_id = (('http://download.jensenlab.org/human_human_textmining_full.tsv', '9606'),
                           ('http://download.jensenlab.org/mouse_mouse_textmining_full.tsv', '10090'),
                           ('http://download.jensenlab.org/rat_rat_textmining_full.tsv', '10116'),
                           ('http://download.jensenlab.org/yeast_yeast_textmining_full.tsv', '4932'))
text_mining_combined_file = os.path.join(data_dir, 'rain_textmining.tsv.gz')
# remove text mining results to make sure they are downloaded before each run as they are continuously updated
if os.path.exists(text_mining_combined_file):
    os.remove(text_mining_combined_file)

string_aliases_100 = 'http://string-db.org/newstring_download/protein.aliases.v10.txt.gz'
string_aliases_91 = 'http://string91.embl.de/newstring_download/protein.aliases.v9.1.txt.gz'
string_aliases_83 = 'http://string83.embl.de/newstring_download/protein.aliases.v8.3.txt.gz'
string_aliases_urls = (string_aliases_100, string_aliases_91, string_aliases_83)

string_species_91 = 'http://string91.embl.de/newstring_download/species.v9.1.txt'
string_species_100 = 'http://string-db.com/newstring_download/species.v10.txt'
string_species_urls = (string_species_100, string_species_91)

gene2ensemble_url = 'ftp://ftp.ncbi.nih.gov/gene/DATA/gene2ensembl.gz'
gene2ensemble_path = os.path.join(data_dir, 'gene2ensembl.gz')

prediction_file_url = 'http://rth.dk/resources/rain/data/raw_predictions'
prediction_file_names = ('starmirdb.tar.gz', 'miRanda_v3.3a.tsv.gz', 'miRDB_v5.0.tsv.gz', 'PITA.tsv.gz', 'RNA22.tsv.gz',
                         'RNAhybrid_seed.tsv.gz', 'targetscan.mammals.tsv.gz')

download_tuples = [
#    (os.path.join(id_dict_dir, 'ncRNAaliasfile.tsv.gz'),
#        'http://rth.dk/resources/rain/pipeline_data/ncRNAaliasfile.tsv.gz'),
    (os.path.join(id_dict_dir, 'ncRNAorthfile.tsv.gz'),
        'http://rth.dk/resources/rain/pipeline_data/ncRNAorthfile.tsv.gz'),
    (os.path.join(data_dir, 'knowledge_channel.tar.gz'),
        'http://rth.dk/resources/rain/pipeline_data/knowledge_channel.tar.gz'),
    (miRTarBase_file, miRTarBase_URL),
    (NPInter_file, NPInter_URL),
    (os.path.join(id_dict_dir, 'NPINTERmapfiles.zip'),
        'http://rth.dk/resources/rain/pipeline_data/NPINTERmapfiles.zip'),
    (os.path.join(id_dict_dir, 'NONCODEID_conversions.tsv'),
        'http://rth.dk/resources/rain/pipeline_data/NONCODEID_conversions.tsv'),
    (mirbase_aliases_txt_path, mirbase_aliases_txt_url),
    (mirbase_mirna_dat_path, mirbase_mirna_dat_url),
    (mir_mature_file_path, mir_mature_file_url),
    (mir_fam_path, mir_fam_url),
    (gene2ensemble_path, gene2ensemble_url),
    (os.path.join(data_dir, 'starbase.tar.gz'),
        'http://rth.dk/resources/rain/pipeline_data/starbase.tar.gz'),
    (os.path.join(data_dir, 'clash.tar.gz'),
        'http://rth.dk/resources/rain/pipeline_data/clash.tar.gz'),
    (os.path.join(data_dir, 'crac.tsv.gz'),
        'http://rth.dk/resources/rain/pipeline_data/crac.tsv.gz'),
    (os.path.join(id_dict_dir, 'ensembl72_protein_lookup.tsv.gz'),
        'http://rth.dk/resources/rain/pipeline_data/ensembl72_protein_lookup.tsv.gz'),
    (os.path.join(id_dict_dir, 'snoRNA_lookup.tsv.gz'),
        'http://rth.dk/resources/rain/pipeline_data/snoRNA_lookup.tsv.gz'),
    (os.path.join(data_dir, 'predictions_precomputed.tsv.gz'),
        'http://rth.dk/resources/rain/pipeline_data/predictions_precomputed.tsv.gz')
]

md5_hash = {
"id_dictionaries/ncRNAorthfile.tsv.gz" : "bc733907d7b26a22505cd67cfb446a41",
"data/knowledge_channel.tar.gz" : "b5251850627c52b268209bce589fce01",
"data/miRTarBase_6.1_MTI.xlsx" : "8235d1db9a9c3a7ff30388a4e61dd1a2",
"data/interaction_NPInter[v2.0].txt.tar.gz" : "8a4e5ad0844b9f4ccd28b91aeb5ba38f",
"id_dictionaries/NPINTERmapfiles.zip" : "741b9cc2433635f84dbfba052e4d3596",
"id_dictionaries/NONCODEID_conversions.tsv" : "c24f88bf24a5ace4c4796cdc519732dc",
"id_dictionaries/miRBase_20_aliases.txt.gz" : "9fcef99ea6896dca4c5c01d2b74504ff",
"id_dictionaries/miRBase_20_miRNA.dat.gz" : "a2f2ec3334638020ddf3cc6c0a1aae71",
"id_dictionaries/mature.fa.gz" : "b1c9f5a61d35da9865e316775f8546ee",
"data/miRBase_20_miFam.dat.gz" : "083e0bd787800834ffb10153aad03de2",
"data/gene2ensembl.gz" : "a17a88a747b5104a4ad91f22971cb565",
"data/starbase.tar.gz" : "d310965db6a79da3272d12717c8dd065",
"data/clash.tar.gz" : "343a304843ba13327f215a4b83658e89",
"data/crac.tsv.gz" : "9b08cd3d044c08ba0c857bdbe4dcc7e3",
"id_dictionaries/ensembl72_protein_lookup.tsv.gz" : "b5139928e2ba489be2f62b5281963ebd",
"id_dictionaries/snoRNA_lookup.tsv.gz" : "579d6f0b6f0790686e63818d6156a769"
}

text_mining_files = []
text_mining_files_tax_ids = []
for text_mining_url, tax in text_mining_urls_tax_id:
    text_mining_file = os.path.join(data_dir, os.path.basename(text_mining_url))
    text_mining_files.append(text_mining_file)
    text_mining_files_tax_ids.append(tax)
    download_tuples.append((text_mining_file, text_mining_url))
for string_alias_url in string_aliases_urls:
    string_aliases_file_path = os.path.join(data_dir, string_alias_url.split('/')[-1])
    download_tuples.append((string_aliases_file_path, string_alias_url))
for string_species_url in string_species_urls:
    string_species_file = os.path.join(data_dir, string_species_url.split('/')[-1])
    download_tuples.append((string_species_file, string_species_url))
for curr_prediction_file_name in prediction_file_names:
    curr_file_path = os.path.join(data_dir, curr_prediction_file_name)
    curr_file_url = os.path.join(prediction_file_url, curr_prediction_file_name)
    download_tuples.append((curr_file_path, curr_file_url))
# RAIN node information files
gene_info_orgs = ('10090', '10116', '4932', '6239', '7227', '7955', '9544', '9606')
for org in gene_info_orgs:
    org_file = org + '_gene_info.html.gz'
    org_path = os.path.join(id_dict_dir, org_file)
    org_url = 'http://rth.dk/resources/rain/pipeline_data/' + org_file
    download_tuples.append((org_path, org_url))

for curr_path, curr_url in download_tuples:
    if curr_path in md5_hash:
        md5 = md5_hash[curr_path]
    else:
        md5 = None
    check_and_download(curr_path, curr_url, force_download=force_dl,
                       skip_unavailable=skip_unavailable_dl, md5sum=md5)

logger.info('Aggregating and compressing species-specific textmining files to: {}'.format(text_mining_combined_file))
for text_mining_file, tax_id in zip(text_mining_files, text_mining_files_tax_ids):
    with open(text_mining_file) as f_in, gzip.open(text_mining_combined_file, 'ab') as f_out:
        for line in f_in:
            f_out.write(tax_id + '\t' + line)
    os.remove(text_mining_file)
