# Read
    library(magrittr)
    library(multicrispr)
    bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
    gr <- bed_to_granges(bedfile, 'mm10')

# Entrez    
    dbname <- 'TxDb.Mmusculus.UCSC.mm10.knownGene'
    db <- utils::getFromNamespace(dbname, dbname)
    annotate_granges(gr, db)    %>% 
    (function(x) x$gene_id)     %>%
    (function(x) x[!is.na(x)])  %>% 
    strsplit(';')               %>% 
    unlist()                    %>% 
    unique()                    %>% 
    write.table('inst/extdata/SRF.entrez', 
                quote = FALSE, row.names = FALSE, col.names = FALSE)

# Ensembl
    db <- EnsDb.Mmusculus.v98()
    annotate_granges(gr, db)    %>% 
    (function(x) x$gene_id)     %>%
    (function(x) x[!is.na(x)])  %>% 
    strsplit(';')               %>% 
    unlist()                    %>% 
    unique()                    %>% 
    write.table('inst/extdata/SRF.ensembl', 
                quote = FALSE, row.names = FALSE, col.names = FALSE)
