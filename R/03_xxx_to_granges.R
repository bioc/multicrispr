
#=============================================================================
# bedfile -> data.table -> GRanges
#=============================================================================

add_seqinfo <- function(gr, bsgenome){
    GenomeInfoDb::seqinfo(gr) <- GenomeInfoDb::seqinfo(bsgenome)
    gr
}


#' Get BSgenome
#' @param gr \code{\link[GenomicRanges]{GRanges-class}}
#' @return BSgenome
#' @examples 
#' bedfile  <- system.file('extdata/SRF.bed', package='multicrispr')
#' gr <- bed_to_granges(bedfile, 'mm10')
#' get_bsgenome(gr)
#' @export
get_bsgenome <- function(gr){
    . <- NULL
    
    assert_is_identical_to_true(is(gr, 'GRanges'))
    genome <- unique(unname(GenomeInfoDb::genome(gr)))
    assert_is_a_string(genome)
    getBSgenome(genome)
}


#' Read bedfile into GRanges
#' 
#' Read bedfile into granges and output graphical and textual summaries
#' 
#' @param bedfile        file path
#' @param genome         genome identifier ('mm10')
#' @param do_order       logical(1)
#' @param plot           logical(1)
#' @param verbose        logical(1)
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @examples
#' bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' bed_to_granges(bedfile, 'mm10')
#' @seealso \code{rtracklayer::import.bed} (documented in 
#' \code{\link[rtracklayer]{BEDFile-class}}), around which this function wraps.
#' @export
bed_to_granges <- function(
    bedfile, 
    genome   = NULL, 
    do_order = TRUE,
    plot     = TRUE, 
    verbose  = TRUE 
){
    . <- NULL
    
    # Assert
    assertive.files::assert_all_are_existing_files(bedfile)
    assertive.types::assert_is_a_bool(plot)
    assertive.types::assert_is_a_bool(verbose)

    # Read
    gr <- rtracklayer::import.bed(bedfile, genome = genome)
    if (verbose) cmessage('\t\t%d ranges on %d chromosomes',
                    length(gr), length(unique(GenomeInfoDb::seqnames(gr))))
    
    # Plot
    title <- paste0(genome, ': ', basename(bedfile))
    if (plot) plot_karyogram(gr, title)
    
    # Order
    if (do_order)
        gr %<>% extract(order(GenomeInfoDb::seqnames(.), start(.)))
    
    # Return
    gr
}

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
#' @param file    Entrez Gene identifier file (one per row)
#' @param geneids Entrez Gene identifier vector
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
    geneids <- utils::read.table(file)[[1]] %>% as.character()
    genes_to_granges(geneids, db, plot = plot)
}


