
#' Annotate GRanges
#' @param granges \code{\link[GenomicRanges]{GRanges-class}}
#' @param db      \code{\link[GenomicFeatures]{TxDb-class}} or 
#'                \code{\link[ensembldb]{EnsDb-class}}
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @seealso \code{\link{genes_to_granges}} (inverse operation)
#' @examples 
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' granges <- read_bed(bedfile, 'mm10', plot = FALSE)
#' db <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
#' annotate_granges(granges, db)
#' @export
annotate_granges <- function(granges, db){
    
    gene_id <- NULL
    
    granges %>% 
    plyranges::join_overlap_left(GenomicFeatures::genes(db)) %>% 
    data.table::as.data.table() %>% 
    extract(!is.na(gene_id) ,  
            gene_id := paste0(gene_id, collapse = ';'), 
            by = c('seqnames', 'start', 'end', 'strand')) %>% 
    unique() %>% 
    data.table::setnames('gene_id', 'is_within') %>% 
    as('GRanges')

}


#' Convert geneids into GRanges
#' @param geneids Entrezg vector
#' @param db    \code{\link[GenomicFeatures]{TxDb-class}} or 
#'                \code{\link[ensembldb]{EnsDb-class}}
#' @param plot TRUE or FALSE
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @seealso \code{\link{annotate_granges}} (inverse operation)
#' @examples
#' # Mouse Entrez
#'     dbname <- 'TxDb.Mmusculus.UCSC.mm10.knownGene'
#'     db <- utils::getFromNamespace(dbname, dbname)
#'     geneids <- names(GenomicFeatures::genes(db))[1:100]
#'     genes_to_granges(geneids, db)
#'
#' # Human Entrez
#'     dbname <- 'TxDb.Hsapiens.UCSC.hg38.knownGene'
#'     db <- utils::getFromNamespace(dbname, dbname)
#'     geneids <- names(GenomicFeatures::genes(db))[1:100]
#'     genes_to_granges(geneids, db)
#'
#' # Mouse Ensembl
#'     db <- EnsDb.Mmusculus.v98()
#'     geneids <- names(GenomicFeatures::genes(db))[1:100]
#'     genes_to_granges(geneids, db)
#'
#' # Human Ensembl
#'     db <- EnsDb.Hsapiens.v98()
#'     geneids <- names(GenomicFeatures::genes(db))[1:100]
#'     genes_to_granges(geneids, db)
#' @export
genes_to_granges <- function(geneids, db, plot = TRUE){
    gr <- GenomicFeatures::genes(db)[geneids]
    if (plot) plot_karyogram(gr)
    gr
}
