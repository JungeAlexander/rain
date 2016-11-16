#!/usr/bin/env python

__author__ = 'Ferhat Alkan'

import gzip
import os
import argparse
from lxml import etree
import logging
logger = logging.getLogger(os.path.basename(__file__))
logging.basicConfig(level=logging.INFO)

parser = argparse.ArgumentParser(description=
                                 """Parses given rna alias file and several BioMart dictionaries in html format
                                  to create a corresponding node_info file with node descriptions and external links.
                                  """)
parser.add_argument('-alias', help='Input alias file to generate node info for', required=True)
parser.add_argument('-info', help='Output node info file. (default=master_files/node_info.tsv)',
                    default='master_files/node_info.tsv')
parser.add_argument('-dicpath', help='Alternative path to download BIOMART dictionaries. (default=id_dictionaries/)',
                    default='id_dictionaries/')


# Reads the alias file and generates the node dictionary together with the set of circBase nodes.
# node dictionary is indexed by organism id
def read_aliases_and_group(alias_file):
    inf = gzip.open(alias_file) if alias_file[-3:] == ".gz" else open(alias_file)
    biomart_dic = {}
    circbase = set()
    for line in inf:
        if line[0] != "#":
            cols = line.rstrip().split('\t')
            if cols[0] not in biomart_dic.keys():
                biomart_dic[cols[0]] = {}
                # id_set[cols[0]] = set()

            if cols[3] == 'circBase' and len(cols[2]) == 16:
                circbase.add((cols[0], cols[1], cols[2]))

            biomart_dic[cols[0]][cols[2]] = cols[1]
            # id_set[cols[0]].add((cols[1],cols[3]))
    inf.close()
    return biomart_dic, circbase


#  Generates the BIOMART node information dic (key=node:value=info_set) for given alias dictionary
#  Informations are tagged by their type and infos dictionary is indexed by organism id and
#  infos[orgid]'s are indexed by node ids
def get_node_info_from_biomart_dics(id_dic):
    infos = {}
    muti_id_rows = 0
    for org in id_dic.keys():
        logger.info('Parsing BIOMART info for '+org)
        infos[org] = {}
        html_file = args.dicpath + '/' + org + '_gene_info.html.gz'
        if not os.path.isfile(html_file):
            logger.warning('no biomart dictionary for ' + org)
            continue
        #  HTML file is parsed below
        inf = gzip.open(html_file) if html_file[-3:] == ".gz" else open(html_file)
        for event, allt in etree.iterparse(inf):
            if allt.tag == 'table':
                for tr in allt:  # each row in table
                    if tr[0].tag == 'th':  # Table headers are parsed to assign tags to different info
                        cols = [th.text.split()[0] for th in tr]
                        cols[0] = 'Ensembl'
                        cols[1] = 'genetype'
                        cols[2] = 'description'
                    elif tr[0].tag == 'td':  # Table rows are parsed and info is retrieved below
                        curids = set()
                        for td in tr:
                            if len(td) > 0 and td[0].text in id_dic[org].keys():  # element with an external link
                                curids.add(id_dic[org][td[0].text])
                            elif len(td) == 0 and td.text in id_dic[org].keys():  # element without a link
                                curids.add(id_dic[org][td.text])

                        if len(curids) == 1:  # If a table row refers to one node id
                            curid = curids.pop()
                            if curid not in infos[org]:
                                infos[org][curid] = set()

                            for i in range(len(tr)):
                                if len(tr[i]) == 0 and tr[i].text is not None:
                                    infos[org][curid].add((cols[i], tr[i].text))
                                if len(tr[i]) > 0:
                                    for trc in tr[i]:
                                        if trc.text is not None and 'href' in trc.attrib.keys():
                                            infos[org][curid].add((cols[i], trc.attrib['href']))
                        elif len(curids) > 1:  # If a table row refers to multiple node ids, report it
                            logger.debug('table row mapped to multiple node IDs' + str(curids))
                            muti_id_rows += 1
        if muti_id_rows > 0:
            logger.info('biomart dictionaries mapped {:d} identifiers to multiple node identifiers.'.format(
                    muti_id_rows))
        inf.close()
    return infos


#  Generates the circBase node information dic (key=node:value=info_set) for given set of circbase nodes
#  Informations are tagged by circBase and infos dictionary is indexed by organism id and
#  infos[orgid]'s are indexed by node ids
def get_circbase_info(circ_set):
    infos = {}
    logger.info('Creating circBase info for all')
    for (org, circ_node, circ_id) in circ_set:
        if org not in infos.keys():
            infos[org] = {}
        if circ_node not in infos[org].keys():
            infos[org][circ_node] = set()
        infos[org][circ_node].add(('circBase', 'http://circbase.org/cgi-bin/singlerecord.cgi?id='+circ_id))
        infos[org][circ_node].add(('genetype','circRNA'))
    return infos


#  Generates the miRBase node information dic (key=node:value=info_set) from given alias file
#  Informations are tagged by miRBase and infos dictionary is indexed by organism id and
#  infos[orgid]'s are indexed by node ids
def get_mirbase_info(alias_file):
    mirnas = set()
    logger.info('Parsing mirBase file')
    mir_file = args.dicpath + '/mature.fa.gz'
    mir_desc = {}
    mirf = gzip.open(mir_file)
    for line in mirf:
        if line[0]=='>':
            cols = line.rstrip().split()
            if cols[1][:5] == 'MIMAT':
                mir_desc[cols[1]] = ' '.join(cols[2:])
    mirf.close()

    inf = gzip.open(alias_file) if alias_file[-3:] == ".gz" else open(alias_file)
    infos = {}
    logger.info('Creating miRBase info for all')
    for line in inf:
        if line[0] != "#":
            org, id, alias, source = line.rstrip().split('\t')
            if org not in infos.keys():
                infos[org] = {}
            if source.lower() == 'mirbase' and len(alias) > 5 and alias[:5] == 'MIMAT':
                if id not in infos[org].keys():
                    infos[org][id] = set()
                infos[org][id].add(('miRBase', 'http://www.mirbase.org/cgi-bin/mature.pl?mature_acc=' + alias))
                infos[org][id].add(('genetype', 'miRNA'))
                if mir_desc.has_key(alias):
                    infos[org][id].add(('description', mir_desc[alias]))
                    mirnas.add(id)
    inf.close()
    # manual curation of this specific example
    curation_ex = 'hsa-mir-196a-1'
    infos['9606'][curation_ex] = set()
    infos['9606'][curation_ex].add(('miRBase', 'http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=MI0000238'))
    infos['9606'][curation_ex].add(('genetype', 'miRNA'))
    infos['9606'][curation_ex].add(('description',  'Homo sapiens miR-196a-1'))
    mirnas.add(curation_ex)
    return infos, mirnas # mirna ids with mirBase descriptions

# Main function
def main(parsed=None):
    #  Get node dics and sets
    biomart, circ = read_aliases_and_group(args.alias)
    # Create info dics
    biomart_info = get_node_info_from_biomart_dics(biomart)
    circ_info = get_circbase_info(circ)
    mir_info, mirnas = get_mirbase_info(args.alias)

    #  Report all the information
    out = open(parsed.info, 'w')
    for org in biomart_info.keys():
        link_replaced_mirs = set()
        logger.info('Writing info for '+org)
        for node, info_set in biomart_info[org].items():
            for info in info_set:
                data = info[1]
                if info[0]=='SGD': # fix SGD link
                    data = info[1].replace('locusS','locus/S')
                if info[0].lower() == 'mirbase' and node in mir_info[org].keys(): # Skip mirbase link if MIMAT present
                    link_replaced_mirs.add(node)
                    continue
                if info[0].lower() == 'description' and node in mirnas:
                    continue
                out.write(org+'\t'+node+'\t'+info[0]+'\t'+data+'\n')

        # Save circ_base info
        if org in circ_info.keys():
            for node, info_set in circ_info[org].items():
                for info in info_set:
                    out.write(org+'\t'+node+'\t'+info[0]+'\t'+info[1]+'\n')

        # Save miRBase info
        if org in mir_info.keys():
            for node, info_set in mir_info[org].items():
                for info in info_set:
                    if info[0] == 'genetype' and node in link_replaced_mirs:
                        continue
                    out.write(org+'\t'+node+'\t'+info[0]+'\t'+info[1]+'\n')

        # manually curated example of tRNA
        curation_ex = 'TRV-AAC1-1'
        curation_ex_org = '9606'
        out.write(curation_ex_org + '\t' + curation_ex + '\t' + 'genetype' + '\t' +
                  'transfer RNA' + '\n')
        out.write(curation_ex_org + '\t' + curation_ex + '\t' + 'description' + '\t' +
                  'transfer RNA valine 24 (anticodon AAC)' + '\n')
        out.write(curation_ex_org + '\t' + curation_ex + '\t' + 'HGNC' + '\t' +
                  'http://www.genenames.org/cgi-bin/gene_symbol_report?hgnc_id=34885' + '\n')
    out.close()

if __name__ == "__main__":
    logger.info("Creating node information file.")
    args = parser.parse_args()
    main(args)
    logger.info('Done.' + os.linesep + os.linesep)