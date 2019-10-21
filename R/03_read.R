
#=============================================================================
# bedfile -> data.table -> GRanges
#=============================================================================

add_seqinfo <- function(gr, bsgenome){
    seqinfo(gr) <- seqinfo(bsgenome)
    gr
}


#' Get BSgenome
#' @param gr \code{\link[GenomicRanges]{GRanges-class}}
#' @return BSgenome
#' @examples 
#' bedfile  <- system.file('extdata/SRF.bed', package='multicrispr')
#' gr <- read_bed(bedfile, 'mm10')
#' get_bsgenome(gr)
#' @export
get_bsgenome <- function(gr){
    . <- NULL
    
    assert_is_identical_to_true(is(gr, 'GRanges'))
    genome <- unique(unname(genome(gr)))
    assert_is_a_string(genome)
    getBSgenome(genome)
}


#' Read bedfile into GRanges
#' 
#' Read bedfile into granges and output graphical and textual summaries
#' 
#' @param bedfile        file path
#' @param genome         character: e.g. 'mm10'
#' @param do_order       logical(1)
#' @param plot           logical(1)
#' @param verbose        logical(1)
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @examples
#' bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' read_bed(bedfile, 'mm10')
#' @seealso \code{rtracklayer::import.bed} (documented in 
#' \code{\link[rtracklayer]{BEDFile-class}}), around which this function wraps.
#' @export
read_bed <- function(
    bedfile, 
    genome, 
    do_order = TRUE,
    plot     = TRUE, 
    verbose  = TRUE 
){
    
    # Assert
    assert_all_are_existing_files(bedfile)
    assert_is_character(genome)
    assert_is_a_bool(plot)
    assert_is_a_bool(verbose)

    # Read
    gr <- rtracklayer::import.bed(bedfile, genome = 'mm10')
    if (verbose) cmessage('\t\t%d ranges on %d chromosomes',
                            length(gr), length(unique(seqnames(gr))))
    
    # Plot
    title <- paste0(genome, ': ', basename(bedfile))
    if (plot) plot_karyogram(gr, title)
    
    # Order
    if (do_order)
        gr %<>% extract(order(seqnames(.), start(.)))
    
    # Return
    gr
}
