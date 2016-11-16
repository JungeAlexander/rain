library(biomaRt)
#listMarts()
ensembl <- useMart('ENSEMBL_MART_ENSEMBL')
#listDatasets(ensembl)
human_ensembl <- useDataset('hsapiens_gene_ensembl', mart = ensembl)
# mouse_ensembl <- useDataset('mmusculus_gene_ensembl', mart=ensembl)
# rat_ensembl <- useDataset('rnorvegicus_gene_ensembl', mart=ensembl)
# yeast_ensembl <- useDataset('scerevisiae_gene_ensembl', mart=ensembl)
# TODO add organism specific DBs: MGI, SGD - or just modify which attributes to export in other organisms?

filters <- listFilters(human_ensembl)

attributes <- listAttributes(human_ensembl)
#paste(attributes[attributes$page == 'feature_page', 'name'], collapse = '\', \'')
human_attributes <-
  c(
    'ensembl_gene_id',
    'ensembl_transcript_id',
    'description',
    'transcript_gencode_basic',
    'external_gene_name',
    'external_gene_source',
    'external_transcript_name',
    'external_transcript_source_name',
    'gene_biotype',
    'transcript_biotype',
    'source',
    'transcript_source',
    'status',
    'transcript_status',
    'embl',
    'ens_hs_transcript',
    'ens_hs_translation',
    'entrezgene',
    'entrezgene_transcript_name',
    'genedb',
    'ottg',
    'ottt',
    'hgnc_id',
    'hgnc_symbol',
    'hgnc_transcript_name',
    'mirbase_accession',
    'mirbase_id',
    'mirbase_transcript_name',
    'refseq_ncrna',
    'refseq_ncrna_predicted',
    'rfam',
    'rfam_transcript_name',
    'rnacentral',
    'ucsc'
  )

# below were selected after consulting: http://www.ensembl.org/Help/Faq?id=468 and running:
# filterOptions("biotype", human_ensembl)
ncrna_biotypes <- c('3prime_overlapping_ncRNA',
                    'antisense',
                    'bidirectional_promoter_lncRNA',
                    'lincRNA',
                    'macro_lncRNA',
                    'miRNA',
                    'misc_RNA',
                    'Mt_rRNA',
                    'Mt_tRNA',
                    'non_coding',
                    'polymorphic_pseudogene',
                    'processed_pseudogene',
                    'processed_transcript',
                    'pseudogene',
                    'ribozyme',
                    'rRNA',
                    'scaRNA',
                    'scRNA',
                    'sense_intronic',
                    'sense_overlapping',
                    'snoRNA',
                    'snRNA',
                    'sRNA',
                    'transcribed_processed_pseudogene',
                    'transcribed_unitary_pseudogene',
                    'transcribed_unprocessed_pseudogene',
                    'unitary_pseudogene',
                    'unprocessed_pseudogene',
                    'vaultRNA')

human <- getBM(attributes = human_attributes, filters='biotype', values=ncrna_biotypes, mart = human_ensembl)
human <- getBM(attributes = c('hgnc_symbol', 'rnacentral', 'ensembl_gene_id'),
               filters='biotype', values=ncrna_biotypes, mart = human_ensembl)
