

#' Add inverse strand for each range
#'
#' @param gr \code{\link[GenomicRanges]{GRanges-class}}
#' @param plot     logical(1)
#' @param verbose  logical(1)
#' @return \code{\link[GenomicRanges]{GRanges-class}}, twice as long as input
#' @examples
#' require(magrittr)
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' gr <- bed_to_granges(bedfile, bsgenome, plot = FALSE)
#' add_inverse_strand(gr)
#' @export
add_inverse_strand <- function(gr, plot = TRUE, verbose = TRUE){
    complements <- invertStrand(gr)
    newranges <- c(gr, complements)
    txt <- sprintf('\t\t%d ranges after adding inverse strands',
                    length(newranges))
    if (plot){
        plot_intervals(
            GRangesList(original = gr, complements = complements),
            title = txt)
    }
    if (verbose) cmessage(txt)
    newranges
}

#=============================================================================
# Find cas9 ranges
#=============================================================================

#' Find cas9 sites in targetranges
#' @param gr        \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome  \code{\link[BSgenome]{BSgenome-class}}
#' @param inclcompl logical(1): include complementary strands in search?
#' @param verbose   logical(1): report verbosely?
#' @return   \code{\link[GenomicRanges]{GRanges-class}}
#' @examples
#' bedfile  <- system.file('extdata/SRF.bed', package='multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' gr <- extend(bed_to_granges(bedfile, bsgenome), plot = FALSE)
#' find_cas9s(gr)
#' @export 
find_cas9s <- function(gr, bsgenome, inclcompl = TRUE, verbose = TRUE){

    # Assert
    assertive.types::assert_is_all_of(gr, 'GRanges')
    assertive.types::assert_is_all_of(bsgenome, 'BSgenome')
    assertive.types::assert_is_a_bool(verbose)
    
    # Add complementary strands
    if (verbose) message('\tFind N{20}NGG cas9seqs')
    if (inclcompl) gr %<>% add_inverse_strand(plot = FALSE, verbose = verbose)
    
    # Comply
    start <- substart <- cas9start <- NULL
    end <- subend <- cas9end <- strand <- seqnames <- NULL
    
    # Find cas9s in targetranges
    targetdt <- data.table::as.data.table(gr)
    targetdt [ , seqs := seqs(gr, bsgenome) ]
    res <- targetdt$seqs %>% stringi::stri_locate_all_regex('[ACGT]{21}GG')
    cextract1 <- function(y) y[, 1] %>% paste0(collapse=';')
    cextract2 <- function(y) y[, 2] %>% paste0(collapse=';')
    targetdt [ , substart := vapply( res, cextract1, character(1)) ]
    targetdt [ , subend   := vapply( res, cextract2, character(1)) ]
    
    # Rm cas9-free targetranges
    idx <- targetdt[, substart == 'NA']
    if (sum(idx)>0){
        if (verbose)  cmessage('\t\tRm %d ranges with no cas9sites', sum(idx)) 
        targetdt %<>% extract(!idx)
    }

    # Transform into cas9ranges
    cas9ranges <- tidyr::separate_rows(targetdt, substart, subend) %>%
        data.table::data.table()                                   %>% 
        extract(, substart := as.numeric(substart))                %>% 
        extract(, subend   := as.numeric(subend))                  %>% 
        extract(, seqs     := substr(seqs, substart, subend))      %>%
        extract( strand=='+', cas9start := start + substart - 1  ) %>% 
        extract( strand=='+', cas9end   := start + subend   - 1  ) %>%
        extract( strand=='-', cas9start := end   - subend   + 1  ) %>%
        extract( strand=='-', cas9end   := end   - substart + 1  ) %>%
        extract(, list( seqnames = seqnames, start = cas9start, 
                        end = cas9end,  strand  = strand,  seqs = seqs)) %>% 
        unique() %>% as('GRanges')
    seqinfo(cas9ranges) <- seqinfo(gr)

    # Return
    if (verbose)   cmessage('\t\t%d cas9 seqs across %d ranges', 
                            length(unique(seqs(cas9ranges))), 
                            length(cas9ranges))
    return(cas9ranges)
}

