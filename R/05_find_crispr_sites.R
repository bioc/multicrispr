


#' Find crispr sites in targetranges
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
#' # Find crispr sites
#'     find_crispr_sites(targets)
#' @export 
find_crispr_sites <- function(targets, plot = TRUE, verbose = TRUE){

    # Assert. Import. Comply
    assertive.types::assert_is_all_of(targets, 'GRanges')
    assertive.sets::assert_is_subset('seq',names(GenomicRanges::mcols(targets)))
    assertive.types::assert_is_character(targets$seq)
    assertive.types::assert_is_a_bool(verbose)
    extract <- magrittr::extract
    start <- substart <- crispr_start <- NULL
    end <- subend <- crispr_end <- strand <- seqnames <- NULL
    
    # Find crispr sites in targetranges
    pattern <- '[ACGT]{21}GG'
    if (verbose) cmessage('\tFind %s crispr sites', pattern)
    targetdt <- data.table::as.data.table(targets)
    res <- targetdt$seq %>% stringi::stri_locate_all_regex(pattern)
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
                data.table::data.table()                                    %>% 
                extract(, substart := as.numeric(substart))                 %>% 
                extract(, subend   := as.numeric(subend))                   %>% 
                extract(, seq      := substr(seq, substart, subend))        %>%
                extract( strand=='+', crispr_start := start + substart - 1) %>% 
                extract( strand=='+', crispr_end   := start + subend   - 1) %>%
                extract( strand=='-', crispr_start := end   - subend   + 1) %>%
                extract( strand=='-', crispr_end   := end   - substart + 1) %>%
                extract(, list( seqnames = seqnames, start = crispr_start, 
                                end = crispr_end,  strand = strand,  seq = seq))
    sites <- GenomicRanges::GRanges(unique(sites_dt), 
                                    seqinfo =  GenomeInfoDb::seqinfo(targets))

    # Plot. Message. Return
    if (plot) plot_karyogram(GenomicRanges::GRangesList(
                                target = targets, crispr_site = sites))
    if (verbose)   cmessage('\t\t%d cas9 seqs across %d ranges', 
                        length(unique(sites$seq)), length(sites))
    return(sites)
}

