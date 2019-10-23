
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
    assert_all_are_existing_files(bedfile)
    assert_is_a_bool(plot)
    assert_is_a_bool(verbose)

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
