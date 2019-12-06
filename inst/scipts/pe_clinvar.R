
# Explore ClinVar dataset
#========================

localfile <- '~/.multicrispr/variant_summary.txt.gz'
if (!file.exists(localfile)){
    utils::download.file(
        url = paste0(   
                'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/', 
                'variant_summary.txt.gz'), 
        destfile = localfile)
}
clinvar <- data.table::fread(localfile)
clinvar %<>% extract(Assembly=='GRCh38')                    # 523 K
clinvar %<>% extract(OriginSimple == 'germline')            # 481 K
clinvar[, sort(table(ClinicalSignificance))]
clinvar %<>% extract(ClinicalSignificance == 'Pathogenic')  #  64 K
clinvar[,  c('Assembly', 
             'ClinicalSignificance', 'ClinSigSimple', 
             'Origin', 'OriginSimple', 
             'PhenotypeIDS', 'VariationID', 'SubmitterCategories', 'OtherIDs', 
             'Guidelines', 'TestedInGTR', 'LastEvaluated', 
             'RCVaccession', 'HGNC_ID', 'RS# (dbSNP)', 
             'nsv/esv (dbVar)', 'ChromosomeAccession') := NULL ]
clinvar[, sort(table(Type))]
clinvar[, sort(table(ReviewStatus))]
clinvar[, length(unique(GeneID))]     # 3 790 genes
clinvar[, .(nalleles = .N), by = 'GeneID'][, .(freq = .N), by = nalleles][order(freq)]
clinvar
