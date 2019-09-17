
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
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' gr <- read_bed(bedfile, bsgenome)
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
#' Reads bedfile into GRanges, adds seqinfo(bsgenome), converts 0-based 
#' into 1-based representation, removes duplicates, and visualizes through 
#' karyoplots.
#' 
#' This function accepts bedfiles with metadata (first 6 columns specify
#' ranges, additional columns specify metadata), which the otherwise similar
#' function \code{\link[rtracklayer]{import.bed} doesn not (which was the 
#' motivation to create this function)
#' 
#' @param bedfile        file path
#' @param bsgenome       BSgenome, e.g. BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' @param zero_based     logical(1): whether bedfile GRanges are 0-based
#' @param rm_duplicates  logical(1)
#' @param plot           logical(1)
#' @param verbose        logical(1)
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @note By convention BED files are 0-based. GRanges are always 1-based.
#'       A good discussion on these two alternative codings is given 
#'       by Obi Griffith on Biostars: https://www.biostars.org/p/84686/
#' @examples
#' bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' (gr <- read_bed(bedfile, bsgenome))
#' @importFrom  data.table  :=
#' @export
read_bed <- function(
    bedfile, 
    bsgenome,
    zero_based    = TRUE,
    rm_duplicates = TRUE, 
    plot          = TRUE,
    verbose       = TRUE 
){
    # Assert
    assert_all_are_existing_files(bedfile)
    assert_is_a_bool(verbose)
    assert_is_a_bool(rm_duplicates)
    assert_is_a_bool(zero_based)

    # Comply
    seqnames <- start <- end <- strand <- .N <- gap <- width <- NULL

    # Read
    if (verbose) cmessage('\tRead %s', bedfile)
    dt <- data.table::fread(bedfile, select = c(seq_len(3), 6),
            col.names = c('seqnames', 'start', 'end', 'strand'))
    data.table::setorderv(dt, c('seqnames', 'start', 'end', 'strand'))
    
    # Transform coordinates: 0-based -> 1-based
    if (zero_based){
        if (verbose)    cmessage('\t\tConvert 0 -> 1-based')
        dt[, start := start + 1]
    }
    
    if (verbose) cmessage('\t\tRanges: %d ranges on %d chromosomes', 
                            nrow(dt), length(unique(dt$seqnames)))
    
    # Drop duplicates
    if (rm_duplicates){
        is_duplicated <- cduplicated(dt)
        if (any(is_duplicated)){
            if (verbose) cmessage('\t\t        %d after removing duplicates')
            dt %<>% extract(!duplicated)
        }
    }
        
    # Turn into GRanges
    gr <-  add_seqinfo(as(dt, 'GRanges'), bsgenome)
    
    # Plot and return
    title <- paste0(providerVersion(bsgenome), ': ', basename(bedfile))
    if (plot) plot_karyogram(gr, title)
    gr
}
