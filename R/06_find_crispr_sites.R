


#' Find crispr sites in targetranges
#' @param gr        \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome  \code{\link[BSgenome]{BSgenome-class}}
#' @param plot      TRUE (default) or FALSE
#' @param verbose   TRUE (default) or FALSE
#' @return   \code{\link[GenomicRanges]{GRanges-class}}
#' @examples
#' # Read bed into granges and extend
#'     require(magrittr)
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'     bedfile  <- system.file('extdata/SRF.bed', package='multicrispr')
#'     gr <- bed_to_granges(bedfile, 'mm10')
#'     gr %<>% extend(plot = FALSE)
#'     
#' # Find crispr sites
#'     find_crispr_sites(gr, bsgenome)
#' @seealso \code{\link{find_prime_sites}} to find prime editing sites
#' @export 
find_crispr_sites <- function(gr, bsgenome, plot = TRUE, verbose = TRUE){

    # Assert. Import. Comply
    assert_is_all_of(gr, 'GRanges')
    assert_is_all_of(bsgenome, 'BSgenome')
    assert_is_a_bool(verbose)
    start <- substart <- crispr_start <- NULL
    end <- subend <- crispr_end <- strand <- seqnames <- NULL
    gr %<>% add_seq(bsgenome)
    
    # Find crispr sites in targetranges
    pattern <- '[ACGT]{21}GG'
    if (verbose) cmessage('\tFind %s crispr sites', pattern)
    targetdt <- as.data.table(gr)
    res <- targetdt$seq %>% stri_locate_all_regex(pattern)
    cextract1 <- function(y) y[, 1] %>% paste0(collapse=';')
    cextract2 <- function(y) y[, 2] %>% paste0(collapse=';')
    targetdt [ , substart := vapply( res, cextract1, character(1)) ]
    targetdt [ , subend   := vapply( res, cextract2, character(1)) ]
    
    # Rm crispr-free targetranges
    idx <- targetdt[, substart == 'NA']
    if (sum(idx)>0){
        if (verbose)  cmessage('\t\tRm %d ranges with no crispr site', sum(idx))
        targetdt %<>% extract(!idx)
    }

    # Transform into crispr ranges
    sites_dt  <-  tidyr::separate_rows(targetdt, substart, subend) %>%
                data.table()                                                %>% 
                extract(, substart := as.numeric(substart))                 %>% 
                extract(, subend   := as.numeric(subend))                   %>% 
                extract(, seq      := substr(seq, substart, subend))        %>%
                extract( strand=='+', crispr_start := start + substart - 1) %>% 
                extract( strand=='+', crispr_end   := start + subend   - 1) %>%
                extract( strand=='-', crispr_start := end   - subend   + 1) %>%
                extract( strand=='-', crispr_end   := end   - substart + 1) %>%
                extract(, list( seqnames = seqnames, start = crispr_start, 
                                end = crispr_end,  strand = strand,  seq = seq))
    sites <- GRanges(unique(sites_dt), seqinfo =  seqinfo(gr))

    # Plot. Message. Return
    if (plot){
        gr$set <- 'target'
        sites$set <- 'crisprsite'
        plot_intervals(c(gr, sites), color_var = 'set')
        sites$set <- NULL
    }
    if (verbose)   cmessage('\t\t%d cas9 seqs across %d ranges', 
                        length(unique(sites$seq)), length(sites))
    return(sites)
}

