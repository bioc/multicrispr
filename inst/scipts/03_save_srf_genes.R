# Read
    library(magrittr)
    library(multicrispr)
    bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
    gr <- bed_to_granges(bedfile, 'mm10')

# Functions
    write <- function(x, file){
        write.table(x, file, quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
    
# Entrez    
    dbname <- 'TxDb.Mmusculus.UCSC.mm10.knownGene'
    db <- utils::getFromNamespace(dbname, dbname)
    entrez  <-  annotate_granges(gr, db) %>% 
                magrittr::extract(!is.na(.$gene_id)) %>%
                data.table::as.data.table() %>% 
                magrittr::extract(, 
                          list(seqnames, start, end, strand, entrez=gene_id))
    
# Ensembl
    db <- AnnotationHub::AnnotationHub()[["AH75036"]] #EnsDb.Mmusculus.v98
    ensembl  <- annotate_granges(gr, db) %>% 
                magrittr::extract(!is.na(.$gene_id)) %>%
                data.table::as.data.table() %>%
                magrittr::extract(, 
                          list(seqnames, start, end, strand, ensembl=gene_id))
    
# Intersect
    both <- merge(entrez, ensembl, by = c('seqnames', 'start', 'end', 'strand'))
    set.seed(3)
    idx <- sample(seq_len(nrow(both)), 10)
    both$entrez[ idx] %>% write('inst/extdata/SRF.entrez')
    both$ensembl[idx] %>% write('inst/extdata/SRF.ensembl')
    
