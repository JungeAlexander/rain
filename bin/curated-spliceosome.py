
import urllib
import os

import stringrnautils

################################################################################
# globals
################################################################################
organisms = ['9606', '10090', '10116', '4932']
UNIPROT_QUERY = "http://www.uniprot.org/uniprot/?query=organism:{0}&columns=id,go-id&format=tab"
# UNIPROT_QUERY = "http://www.uniprot.org/uniprot/?query=organism:{0}&columns=id,go-id,go(cellular component)&format=tab"

DATA_DIR = os.path.join(os.path.dirname(__file__), '../data')
goterm_to_urna = {

    ########################################
    # Biological Process (3'-]end processing
    'GO:0034473': ['U1'],
    'GO:0034474': ['U2'],
    'GO:0034475': ['U4'],
    'GO:0034476': ['U5'],
    'GO:0034477': ['U6'],

    ########################################
    # molecular function (binding)

    # uncommented because it does not discriminate with the U rna and the protein
    #'GO:1990446' : 'U1 snRNP',
    #'GO:1990447' : 'U2 snRNP',

    'GO:0030619': ['U1'],
    'GO:0030620': ['U2'],
    'GO:0030621': ['U4'],
    'GO:0030623': ['U5'],
    'GO:0017070': ['U6'],
    'GO:0071209': ['U7'],
    'GO:0030622': ['U4atac'],
    'GO:0030624': ['U6atac'],
    'GO:0030625': ['U11'],
    'GO:0030626': ['U12'],

    ########################################
    # Cellular Component (in complex)
    'GO:0005685': ['U1'],
    'GO:0005686': ['U2'],

    'GO:0005687': ['U4'],
    'GO:0005682': ['U5'],
    'GO:0005688': ['U6'],
    'GO:0005683': ['U7'],

    'GO:0005692': ['U11'],
    'GO:0005693': ['U12'],
    'GO:0005690': ['U4atac'],
    'GO:0005691': ['U6atac'],


    'GO:0071001': ['U4', 'U6'],
    'GO:0034693': ['U11', 'U12'],
    'GO:0071002': ['U4atac', 'U6atac'],
    ########################################
    # 'GO:0071024' : 'SL snRNP',
}



################################################################################
# Functions
################################################################################
def get_u_snrna_rna_rna_interactions():
    interactions = []

    # Core component of the spliceosomal U1, U2, U4 and U5 small nuclear ribonucleoproteins (snRNPs),
    # the building blocks of the spliceosome
    core = ('U1', 'U2', 'U4', 'U5')
    for i, rna1 in enumerate(core, 1):
        for rna2 in core[i:]:
            # a better refrence may be: Maniatis, T. and Reed, R. (1987). Nature 325, 673-678.
            # but I do not have journal acces at my new job and have therefore not read it -jcr
            url = "http://bioscience.jbpub.com/cells/MBIO5245.aspx"
            interactions.append((rna1, rna2, url))

    # splisosomal tri complex GO:0097526
    url = "http://www.ebi.ac.uk/QuickGO/GTerm?id={}"
    interactions.append(('U4', 'U5', url.format('GO:0046540')))
    interactions.append(('U4', 'U6', ','.join((url.format('GO:0046540'), url.format('GO:0071001')))))
    interactions.append(('U5', 'U6', url.format('GO:0046540')))
    interactions.append(('U11', 'U12', url.format('GO:0034693')))
    interactions.append(('U4atac', 'U6atac', url.format('GO:0071002')))
    return interactions


def download():
    for organism in organisms:
        url = UNIPROT_QUERY.format(organism).replace(' ', '%20')
        file_name = os.path.join(DATA_DIR, '{}_celluar_component.tsv'.format(organism))
        print (url)
        print (file_name)
        urllib.urlretrieve(url, file_name)


def run():
    # TODO: UNCOMMENT BEFORE YOU COMMIT!!!
    # download()

    master_file = open ('master_files/database_spliceosome.tsv', 'w')

    # TODO: map RNAs as well!!!
    for organism in organisms:
        print (' - generating interactions for {}'.format(organism))
        interactions = get_u_snrna_rna_rna_interactions()
        for ent1, ent2, url in interactions:
            _str = "\t".join((organism, ent1, ent2, "0", "DATABASE", "0.9", url, ''))
            master_file.write("{}\n".format(_str))

        uniprot_to_ensp = stringrnautils.get_alias_to_string_mapper(organism, '', '', 10)[organism]
        with open('data/{}_celluar_component.tsv'.format(organism)) as go_terms_file:
            go_terms_file.readline()  # skip header
            for line in go_terms_file:
                try:
                    uniprot_acc, goterms = line.rstrip().split('\t')
                except ValueError:
                    # there are no go terms for this protein
                    continue
                url = 'http://www.uniprot.org/uniprot/{}'.format(uniprot_acc)
                for go_term in goterms.split('; '):
                    ensp = uniprot_to_ensp.get(uniprot_acc, None)
                    if ensp:
                        for u_rna in goterm_to_urna.get(go_term, []):
                            _line = "\t".join((organism, u_rna, ensp, "0", "DATABASE", "0.9", url, ''))
                            master_file.write("{}\n".format(_line))


if __name__ == '__main__':
    run()


