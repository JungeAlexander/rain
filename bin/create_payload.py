import collections
import os
import re
import argparse
import logging
logger = logging.getLogger(os.path.basename(__file__))
logging.basicConfig(level=logging.INFO)


from stringrnautils import get_string_to_alias_mapper


def is_valid_file(path):
    if not os.path.isfile(path):
        parser.error("Given node info file {} does not exist.".format(path))
    return path


parser = argparse.ArgumentParser(description=
                                 """Create payload files for STRING.""")
parser.add_argument('-data_path',
                    help='path to data',
                    default = 'data/')
parser.add_argument('-nodes_file',
                    help='path to file with comments for nodes',
                    default = None)
parser.add_argument('-payload_server',
                    help='server where the payload files will be stored',
                    default = 'http://download.jensenlab.org/STRING-RNA/')
parser.add_argument('-project_name',
                    help='Name given to the project',
                    default = 'string-rna payload')
parser.add_argument('-ncrna_info_file', type=is_valid_file,
                    help='tab-delimited file with four columns: '
                         'taxonomy ID, ncRNA identifier, information type (e.g., genetype, description or URL), '
                         'information')
parser.add_argument('-tax_identifiers', type=str, nargs='+',
                    help='taxonomy identifiers of organisms to be used')
args = parser.parse_args()


MASTER_FILE = args.data_path+"all.tsv"

def create_pyload_files(data_path, project_name):
        
    payload_files = [data_path+"edges1"+".tsv",data_path+"edges2"+".tsv",data_path+"nodes"+".tsv",data_path+project_name+".json"]
    
    return payload_files

#JSON format
#{
#"extension_nodes_file" : "http://jensenlab.org/ljj/payloadtest/xnodes1.tsv",
#"edges_file"           : "http://jensenlab.org/ljj/payloadtest/edges1.tsv",
#"extension_edges_file" : "http://jensenlab.org/ljj/payloadtest/xedges1.tsv",
#"name"                 : "Payload test"
#}
def create_json_file(nodes_edges_files,server_url,project_name):
    import json
    
    json_file = nodes_edges_files.pop()
    
    f = open(json_file,'w')
    f.write(json.dumps({'extension_nodes_file': server_url+nodes_edges_files.pop(), 'extension_edges_file':server_url+nodes_edges_files.pop(),
                'edges_file': server_url+nodes_edges_files.pop(),'name': project_name},
        sort_keys=True, indent=4, separators=(',', ': ')))
    f.close()

#Read new nodes extra information
#This information is shown in the iframe when clicking in a new node
#The input file should have the format:
#   identifier    comment    linkout (if necessary)
def parse_nodes_extra_information(nodes_file):
    nodes_extra_info = collections.defaultdict(dict)
    if nodes_file is not None:
        with open(nodes_file,'r') as nf:
            for line in nf:
                data = line.strip().split("\t")
                nodes_extra_info[data[0]] = data[1]+"\t"+data[2]
        nf.close()
        
    return nodes_extra_info


def parse_ncrna_info_file(ncrna_info_file_path):
    """
    :param ncrna_info_file_path: path to tab-delimited file with four columns:
           taxonomy ID, ncRNA identifier, information type (e.g., genetype, description or URL), information
    :return: a dict mapping ncRNA identifiers to comment and linkout information to be displayed in the network,
             the two columns separated by a tab character
    """
    reference_db_names = {'circBase', 'Ensembl', 'FlyBase', 'HGNC', 'MGI', 'miRBase', 'SGD', 'WormBase', 'ZFIN'}

    ncrna_identifiers = set()
    ncrna_to_genetype = collections.defaultdict(set)
    ncrna_to_description = collections.defaultdict(set)
    ncrna_to_source_to_url = collections.defaultdict(dict)

    with open(ncrna_info_file_path, 'r') as ncrna_info_file_handle:
        for line in ncrna_info_file_handle:
            tax, ncrna, info_type, info = line.rstrip('\n\r').split('\t')
            ncrna_identifiers.add(ncrna)
            if info_type == 'genetype':
                ncrna_to_genetype[ncrna].add(info)
            elif info_type == 'description':
                ncrna_to_description[ncrna].add(info)
            elif info_type in reference_db_names:
                if info_type in ncrna_to_source_to_url[ncrna]:
                    # only the first link per reference db is kept
                    continue
                else:
                    ncrna_to_source_to_url[ncrna][info_type] = info
            else:
                raise ValueError('Error when parsing file {}. Unrecognized type of information: {}'.format(
                    ncrna_info_file_path, info_type))
    logger.info('Found information for {:d} ncRNAs.'.format(len(ncrna_identifiers)))

    ncrna_to_columns = {}
    for ncrna in ncrna_identifiers:
        genetype = ''
        if ncrna in ncrna_to_genetype:
            genetype = 'genetype: ' + '; '.join(sorted(list(ncrna_to_genetype[ncrna])))
        description = ''
        if ncrna in ncrna_to_description:
            description = 'description: ' + '; '.join(sorted(list(ncrna_to_description[ncrna])))
        comment_column = description + '<br>' + genetype

        linkout_column = ''
        sep = ''
        if ncrna in ncrna_to_source_to_url:
            for ref_db, url in ncrna_to_source_to_url[ncrna].iteritems():
                linkout_column += sep + ref_db + '|' + url
                sep = '|'

        ncrna_to_columns[ncrna] = comment_column + '\t' + linkout_column
    return ncrna_to_columns

#Master file format
#Organism        Id1     Id2     Directed        Evidence        Score   Source  URL
def parse_ncRNA_master_file(string_data, payload_files,payload_server,project_name, nodes_extra_info):
    import csv
    
    new_edges1 = open(payload_files[0],'w')
    new_edges2 = open(payload_files[1],'w')
    new_nodes = open(payload_files[2],'w')

    organisms_not_in_list = set()
    info_nodes = collections.defaultdict(set)
    saved_nodes = collections.defaultdict(set)
    for line in open(MASTER_FILE, 'r'):
        if(re.match('^\s*#',line)):
            continue
        row = line.rstrip('\r\n').split('\t')

        try:
            organism, node1, node2, directed, evidence, score, source, linkout, comment = row
        except ValueError as e:
            logger.critical('Error (see below) encountered when unpacking row in aggregate master file. '
                             'Corresponding row was:\n%s' % row)
            raise e
        
        if organism not in taxonomy_ids:
            organisms_not_in_list.add(organism)
            continue

        # use source as comment for edge for now as comment field not populated
        comment = source
        
        node1_extra = "\t"
        node2_extra = "\t"
        if node1 in nodes_extra_info.keys():
            node1_extra = nodes_extra_info[node1]
            info_nodes[organism].add(node1)
        if node2 in nodes_extra_info.keys():
            node2_extra = nodes_extra_info[node2]
            info_nodes[organism].add(node2)

        if(string_data[organism].has_key(node1) and string_data[organism].has_key(node2)):
            line = organism + "\t" + organism+'.'+node1+"\t"+organism+'.'+node2+"\t"+evidence+"\t"+score+"\t"+comment+"\t"+linkout+"\n"
            new_edges1.write(line)
        elif(not string_data[organism].has_key(node1) and string_data[organism].has_key(node2)):
            line = organism + "\t" + organism+'.'+node2+"\t"+node1+"\t"+evidence+"\t"+score+"\t"+comment+"\t"+linkout+"\n"
            new_edges2.write(line)
            if not node1 in saved_nodes[organism]:
                new_nodes.write(organism + "\t" + node1+"\t\t"+node1_extra+"\n")
                saved_nodes[organism].add(node1)
        elif(string_data[organism].has_key(node1) and not string_data[organism].has_key(node2)):
            line = organism + "\t" + organism+'.'+node1+"\t"+node2+"\t"+evidence+"\t"+score+"\t"+comment+"\t"+linkout+"\n"
            new_edges2.write(line)
            if not node2 in saved_nodes[organism]:
                new_nodes.write(organism + "\t" + node2+"\t\t"+node2_extra+"\n")
                saved_nodes[organism].add(node2)
        else:
            line = organism + "\t" + node1+"\t"+node2+"\t"+evidence+"\t"+score+"\t"+comment+"\t"+linkout+"\n"
            new_edges2.write(line)
            if not node1 in saved_nodes[organism]:
                new_nodes.write(organism + "\t" + node1+"\t\t"+node1_extra+"\n")
            if not node2 in saved_nodes[organism]:
                new_nodes.write(organism + "\t" + node2+"\t\t"+node2_extra+"\n")
            saved_nodes[organism].add(node1)
            saved_nodes[organism].add(node2)
    new_edges1.close()
    new_edges2.close()
    new_nodes.close()
    logger.info('Removed interactions for following taxonomy IDs because they were not in'
                'species list: ' + ', '.join(sorted(organisms_not_in_list)))

    node_count = sum([len(v) for k, v in saved_nodes.iteritems()])
    logger.info('Found a total of {:d} new nodes.'.format(node_count))

    info_node_count = sum([len(v) for k, v in info_nodes.iteritems()])
    logger.info('Found node information for {:d} nodes.'.format(info_node_count))

    create_json_file(payload_files,payload_server,project_name)

############################################################################################
if __name__ == '__main__':
    logger.info("Creating payload files.")
    taxonomy_ids = args.tax_identifiers
    ncrna_info_file = args.ncrna_info_file
    logger.info('Retaining following taxonomy identifiers: ' + ' '.join(taxonomy_ids))

    string_data = get_string_to_alias_mapper(taxonomy_ids,'','',10,'all',True)
    payload_files = create_pyload_files(args.data_path,args.project_name)

    if ncrna_info_file:
        ncrna_info_map = parse_ncrna_info_file(ncrna_info_file)
    else:
        ncrna_info_map = {}

    parse_ncRNA_master_file(string_data, payload_files, args.payload_server, args.project_name, ncrna_info_map)
    logger.info('Done.' + os.linesep + os.linesep)
