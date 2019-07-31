
#' Convert range to sequence
#' @param chr     character vector: c('chr1', 'chr1', ...)
#' @param start   numeric vector: start positions
#' @param end     numeric vector: end positions
#' @param strand  character vector: +/- values
#' @param bsgenome  BSgenome object
#' @return character vector: sequences
#' @examples
#' range2seq(
#'     chr      = c('chr2', 'chr1', 'chr17'),
#'     start    = c(83972889, 85603908, 70963622),
#'     end      = c(83972905, 85603924, 70963638),
#'     strand   = c('+', '-', '-'),
#'     bsgenome = BSgenome.Mmusculus.UCSC.mm10::Mmusculus)
#' @export
range2seq <- function(chr, start, end, strand, bsgenome){
    
    # Assert
    assertive.types::assert_is_character(chr)
    assertive.types::assert_is_numeric(start)
    assertive.types::assert_is_numeric(end)
    assertive.types::assert_is_character(strand)
    assertive.sets::assert_is_subset(unique(strand), c('+', '-'))
    tmp <- Reduce( assertive.properties::assert_are_same_length, 
                   list(chr, start, strand))
    assertive.base::assert_is_identical_to_true(methods::is(bsgenome, 'BSgenome'))
    
    # Return
    . <- NULL
    GenomicRanges::GRanges(
        seqnames = chr,
        ranges   = IRanges::IRanges(start = start, end = end),
        strand   = strand,
        seqinfo  = GenomeInfoDb::seqinfo(bsgenome)) %>%
    BSgenome::getSeq(bsgenome, .) %>%
    as.character()
}


#' Find cas9 sites in given ranges
#' @param ranges data.table(chr, start, end, strand)
#' @param bsgenome BSgenome object (eg BSgenome.Mmusculus.UCSC.mm10::Mmusculus)
#' @param verbose logical(1)
#' @return data.table(chr, start, end, strand, cas9seq, cas9start, cas9end). \cr
#'         One row per cas9 site
#' @examples
#' require(magrittr)
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' bedfile  <- system.file('extdata/SRF_sites.bed', package='crisprapex')
#' read_bed(bedfile) %>% slop_fourways(bsgenome) %>% find_cas9s(bsgenome)
#' @export 
find_cas9s <- function(ranges, bsgenome, verbose = TRUE){
    
    # Comply
    chr <- start <- end <- strand <- NULL
    substart <- subend <- cas9start <- cas9end <- cas9seq <- NULL
    
    # Find cas9sites
    if (verbose) message('\tFind N{20}NGG cas9 sites in provided ranges')
    ranges [ , seq := range2seq(chr, start, end, strand, bsgenome) ] 
    res <- ranges$seq %>% stringi::stri_locate_all_regex('[ACGT]{21}GG')
    ranges [ , substart := vapply(  res, 
                                    function(y) y[, 1] %>% paste0(collapse=';'),
                                    character(1)) ]
    ranges [ , subend   := vapply(  res,
                                    function(y) y[, 2] %>% paste0(collapse=';'),
                                    character(1)) ]
    
    # Rm cas9-free ranges
    idx <- ranges[, substart != 'NA']
    if (sum(idx)>0){
        if (verbose)  cmessage('\t\tRm %d ranges with no cas9sites', sum(idx)) 
        ranges %<>% extract(idx)
    }

    # Calculate cas9 ranges
    cas9dt <-   ranges %>% 
                tidyr::separate_rows(substart, subend) %>%
                data.table::data.table() %>% 
                extract(, substart := as.numeric(substart)) %>% 
                extract(, subend   := as.numeric(subend))
    
    cas9dt [ , cas9seq   := substr(seq, substart, subend)    ]
    cas9dt [ , seq := NULL ]
    cas9dt [ strand=='+', cas9start := start + substart - 1  ]
    cas9dt [ strand=='+', cas9end   := start + subend   - 1  ]
    cas9dt [ strand=='-', cas9start := end   - subend   + 1  ]
    cas9dt [ strand=='-', cas9end   := end   - substart + 1  ]
    cas9dt [ , c('substart', 'subend') := NULL]
    
    # Return
    if (verbose)   cmessage('\t\t%d cas9 seqs across %d ranges', 
                            cas9dt[, length(unique(cas9seq))], nrow(cas9dt))
    return(cas9dt)
}


#' Rm center cas9s
#' @param flank_cas9s   data.table(cas9seq, ...)
#' @param center_cas9s  data.table(cas9seq, ...)
#' @param verbose logical(1)
#' @examples
#' require(magrittr)
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' bedfile <- system.file('extdata/SRF_sites.bed', package='crisprapex')
#' tbranges <- read_bed(bedfile)
#' flank9s  <- tbranges %>% flank_fourways(bsgenome) %>% find_cas9s(bsgenome)
#' center9s <- tbranges %>% slop_fourways(bsgenome)  %>% find_cas9s(bsgenome)
#' flank9s %>% rm_center_cas9s(center9s)
#'
#' @return subset of flank_cas9s
#' @export
rm_center_cas9s <- function(flank_cas9s, center_cas9s, verbose = TRUE){
    
    # Assert
    assertive.types::assert_is_data.table(flank_cas9s)
    assertive.types::assert_is_data.table(center_cas9s)
    assertive.types::assert_is_a_bool(verbose)
    
    # Remove
    region <- ncenter <- NULL
    cas9dt <- rbind(cbind(flank_cas9s, region = 'target'), 
                    cbind(center_cas9s,  region = 'taboo'))
    if (verbose)   cmessage(
                    '\t\t%d cas9 seqs across %d ranges', 
                    length(unique(cas9dt$cas9seq)), nrow(cas9dt))
    cas9dt [ , ncenter := sum(region == 'center'), by = 'cas9seq' ]
    cas9dt %<>% extract(ncenter == 0)
    cas9dt [ , c('ncenter', 'region') := NULL ]
    
    # Return
    if (verbose)   cmessage(
                    '\t\t%d cas9 seqs across %d ranges', 
                    length(unique(cas9dt$cas9seq)), nrow(cas9dt))
    cas9dt
}
