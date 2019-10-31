

#=============================================================================
# Find cas9 ranges
#=============================================================================

#' Find cas9 sites in targetranges
#' @param targets   \code{\link[GenomicRanges]{GRanges-class}}
#' @param plot      TRUE (default) or FALSE
#' @param verbose   TRUE (default) or FALSE
#' @return   \code{\link[GenomicRanges]{GRanges-class}}
#' @examples
#' # Read bed into granges, extend, add seqs
#'     require(magrittr)
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'     bedfile  <- system.file('extdata/SRF.bed', package='multicrispr')
#'     targets <- bed_to_granges(bedfile, 'mm10')
#'     targets %<>% extend(plot = FALSE)
#'     targets %<>% add_seq(bsgenome)
#'     
#' # Find cas9 seqs
#'     find_cas9s(targets)
#' @export 
find_cas9s <- function(targets, plot = TRUE, verbose = TRUE){

    # Assert. Comply
    assertive.types::assert_is_all_of(targets, 'GRanges')
    assertive.sets::assert_is_subset('seq', names(mcols(targets)))
    assertive.types::assert_is_character(targets$seq)
    assertive.types::assert_is_a_bool(verbose)
    start <- substart <- cas9start <- NULL
    end <- subend <- cas9end <- strand <- seqnames <- NULL
    
    # Find cas9s in targetranges
    if (verbose) message('\tFind N{20}NGG cas9seqs')
    targetdt <- data.table::as.data.table(targets)
    res <- targetdt$seq %>% stringi::stri_locate_all_regex('[ACGT]{21}GG')
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
    cas9s <- tidyr::separate_rows(targetdt, substart, subend) %>%
        data.table::data.table()                                   %>% 
        extract(, substart := as.numeric(substart))                %>% 
        extract(, subend   := as.numeric(subend))                  %>% 
        extract(, seq      := substr(seq, substart, subend))      %>%
        extract( strand=='+', cas9start := start + substart - 1  ) %>% 
        extract( strand=='+', cas9end   := start + subend   - 1  ) %>%
        extract( strand=='-', cas9start := end   - subend   + 1  ) %>%
        extract( strand=='-', cas9end   := end   - substart + 1  ) %>%
        extract(, list( seqnames = seqnames, start = cas9start, 
                        end = cas9end,  strand  = strand,  seq = seq)) %>% 
        unique() %>% as('GRanges')
    seqinfo(cas9s) <- seqinfo(targets)

    # Plot/Message
    if (plot) plot_karyogram(GenomicRanges::GRangesList(
                                target = targets, cas9site = cas9s))
    if (verbose)   cmessage('\t\t%d cas9 seqs across %d ranges', 
                        length(unique(cas9s$seq)), length(cas9s))
    
    # Return
    return(cas9s)
}

