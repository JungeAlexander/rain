
# core lib
import urllib
import collections
import os.path

# 3rd party
import pandas as pd

# local
import stringrnautils


# Y2 Hybrid (counts as 1 extra pmid)
TABEL_S3_URL = "http://www.cell.com/cms/attachment/2007964525/2030738680/mmc4.xls"
TABEL_S3_SHEET = "TableS3 - Y2H interactions"
TABEL_S3_LOCAL = "data/hegel_s3.xls"

# Litterature Cureted (pmids scores)
TABEL_S4_URL = "http://www.cell.com/cms/attachment/2007964525/2030738681/mmc5.xls"
TABEL_S4_SHEET = "TableS4-literature interactions"
TABEL_S4_LOCAL = "data/hegel_s4.xls"


# luciferase confirmed (0.9 database)
TABEL_S6_URL = "http://www.cell.com/cms/attachment/2007964525/2030738679/mmc6.xls"
TABEL_S6_SHEET = "TableS6 - Co-IP data"
TABEL_S6_LOCAL = "data/hegel_s6.xls"

def run():
    protein_mapper = stringrnautils.get_alias_to_string_mapper("9606", 'ENSP', "")["9606"]


    # parse S3 and S4, use S4 as simply yet another PID supporting them (22365833)


    ############################################################
    # parse S3 and S4 into pmids and benchmark
    ############################################################

    pmid_evicence = collections.defaultdict(list)  # s3 and s4
    experiment_evidence = collections.defaultdict(list)  # s6 counts as 0.9
    hegel_pmid = 22365833

    if not os.path.exists(TABEL_S3_LOCAL):
        urllib.urlretrieve(TABEL_S3_URL, TABEL_S3_LOCAL)
    s3_data = pd.read_excel(TABEL_S3_LOCAL, TABEL_S3_SHEET)

    for i in range(s3_data.shape[0]):
        gene_a = s3_data["SymbolA"][i]
        gene_b = s3_data["SymbolB"][i]

        try:
            prot_a = protein_mapper[gene_a]
            prot_b = protein_mapper[gene_b]

            key = tuple(sorted((prot_a, prot_b)))
            pmid_evicence[key].append(hegel_pmid)
        except KeyError:
            print "the interaction between {0} and {1} could not be mapped to string ids".format(gene_a, gene_b)

    if not os.path.exists(TABEL_S4_LOCAL):
        urllib.urlretrieve(TABEL_S4_LOCAL, TABEL_S4_SHEET)
    s4_data = pd.read_excel(TABEL_S4_LOCAL, TABEL_S4_SHEET)
    for i in range(s4_data.shape[0]):
        gene_a = s4_data["SymbolA"][i]
        gene_b = s4_data["SymbolB"][i]

        try:
            prot_a = protein_mapper[gene_a]
            prot_b = protein_mapper[gene_b]

            key = tuple(sorted((prot_a, prot_b)))
            for j in range(1, int(s4_data["#PMID"][i])):
                pmid_evicence[key].append(int(s4_data[str(j)][i]))

        except KeyError:
            print "the interaction between {0} and {1} could not be mapped to string ids".format(gene_a, gene_b)


    # TODO alex or garde: hook this into the combine-miRTarBase-NPinter script, as thise interactions are proteins and
    # therefore have to be scored using the bins from there (as we have no protein positive set
    for (prot_1, prot_b), pmids in pmid_evicence.items():
        print '\t'.join(("9606", prot_a, prot_b, "0", "Experiment", str(len(pmids)), "Litterature", "", ""))

    ############################################################
    # parse S6 into experiments
    ############################################################
    if not os.path.exists(TABEL_S6_LOCAL):
        urllib.urlretrieve(TABEL_S6_URL, TABEL_S6_LOCAL)
    s6_data = pd.read_excel(TABEL_S6_LOCAL, TABEL_S6_SHEET)

    for i in range(s6_data.shape[0]):
        gene_a = s6_data["FireSymbol"][i]
        gene_b = s6_data["PASymbol"][i]

        try:
            prot_a = protein_mapper[gene_a]
            prot_b = protein_mapper[gene_b]

            key = tuple(sorted((prot_a, prot_b)))
            for j in range(1, int(s4_data["#PMID"][i])):
                pmid_evicence[key].append(int(s4_data[str(j)][i]))

        except KeyError:
            print "the interaction between {0} and {1} could not be mapped to string ids".format(gene_a, gene_b)

        # TODO: append this to the experiments file in "combine_experiments"
        # these are proteins and therefore cannot be benchmarked against our gold-standard
        print '\t'.join(("9606", prot_a, prot_b, "0", "Experiment", "0.9", "Luciferase",
                         "https://www.ncbi.nlm.nih.gov/pubmed/22365833", ""))


if __name__ == '__main__':
    run()