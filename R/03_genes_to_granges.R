
#' Annotate GRanges
#' @param granges \code{\link[GenomicRanges]{GRanges-class}}
#' @param db      \code{\link[GenomicFeatures]{TxDb-class}} or 
#'                \code{\link[ensembldb]{EnsDb-class}}
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @seealso \code{\link{genes_to_granges}} (inverse operation)
#' @examples
#' # Read 
#'     bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'     granges <- bed_to_granges(bedfile, 'mm10', plot = FALSE)
#'     
#' # Annotate Entrez 
#'     dbname <- 'TxDb.Mmusculus.UCSC.mm10.knownGene'
#'     db <- utils::getFromNamespace(dbname, dbname)
#'     annotate_granges(granges, db)
#' 
#' # Annotate Ensembl
#'     db <- EnsDb.Mmusculus.v98()
#'     annotate_granges(granges, db)
#' @export
annotate_granges <- function(granges, db){

    gene_id <- NULL
    if (seqlevelsStyle(granges) != seqlevelsStyle(db)){
        message("Setting seqlevelsStyle(db) <- seqlevelsStyle(granges)")
        seqlevelsStyle(db) <- seqlevelsStyle(granges)
    }
    
    GR <- as(data.table::as.data.table(granges), 'GRanges')
    DR <- GenomicFeatures::genes(db)                                     %>%
            data.table::as.data.table()                                  %>%
            extract(seqlevelsInUse(granges), on = 'seqnames')            %>%
            extract(, c('seqnames', 'start', 'end', 'strand', 'gene_id'),
                with = FALSE)                                            %>%
            as('GRanges')
    
    GR %<>% plyranges::join_overlap_left(DR)                        %>%
            data.table::as.data.table()                             %>% 
            extract(!is.na(gene_id) ,  
                    gene_id := paste0(gene_id, collapse = ';'), 
                    by = c('seqnames', 'start', 'end', 'strand'))   %>%
            unique()                                                %>% 
            as('GRanges')
    GR
}


#' Convert geneids into GRanges
#' @param geneids Entrezg vector
#' @param db    \code{\link[GenomicFeatures]{TxDb-class}} or 
#'                \code{\link[ensembldb]{EnsDb-class}}
#' @param plot TRUE or FALSE
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @seealso \code{\link{annotate_granges}} (inverse operation)
#' @examples
#'  # Entrez
#'      dbname  <-'TxDb.Mmusculus.UCSC.mm10.knownGene'
#'      db      <- utils::getFromNamespace(dbname, dbname)
#'      genefile <- system.file('extdata/SRF.entrez', package='multicrispr')
#'      genefile_to_granges(genefile, db)
#'      geneids <- as.character(read.table(genefile)[[1]])
#'      genes_to_granges(geneids, db)
#'
#' # Ensembl
#'      db <- EnsDb.Mmusculus.v98()
#'      genefile <- system.file('extdata/SRF.ensembl', package='multicrispr')
#'      genefile_to_granges(genefile, db)
#'      geneids <- as.character(read.table(genefile)[[1]])
#'      genes_to_granges(geneids, db)
#'      
#' @export
genes_to_granges <- function(geneids, db, plot = TRUE){
    gr <- GenomicFeatures::genes(db)[geneids]
    if (plot) plot_karyogram(gr)
    gr
}


#' @rdname genes_to_granges
#' @export
genefile_to_granges <- function(file, db, plot = TRUE){
    assertive.files::assert_all_are_existing_files(file)
    geneids <- read.table(file)[[1]] %>% as.character()
    genes_to_granges(geneids, db, plot = plot)
}

