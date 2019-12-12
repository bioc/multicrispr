
cextract1 <- function(y) y[, 1] %>% paste0(collapse=';')
cextract2 <- function(y) y[, 2] %>% paste0(collapse=';')
vextract1 <- function(z) vapply(z, cextract1, character(1))
vextract2 <- function(z) vapply(z, cextract2, character(1))

copy <- function(
    gr, 
    start  = GenomicRanges::start(gr), 
    end    = GenomicRanges::end(gr),
    strand = GenomicRanges::strand(gr), 
    ...
){
    newgr <- gr
    GenomicRanges::start(gr) <- start
    GenomicRanges::end(gr)   <- end
    strand(gr) <- strand
    
    dots_list <- list(...)
    for (mvar in names(dots_list)) mcols(gr)[[mvar]] <- dots_list[[mvar]]
    gr
}

#' Extend ranges to find GG
#' 
#' Extend ranges to find prime editing GG
#' 
#' Extends each target range to the area in which to search for a prime editing 
#' GG duplet, as shown in the sketch below.
#' 
#'                    ===============>
#'                              ----GG--------->
#'                ----GG--------->
#'                              **
#'                <---------GG---
#'                              <---------GG----
#'                          <===============    
#'         
#' @param gr   target \code{\link[GenomicRanges]{GRanges-class}}
#' @param nrt  n RT nucleotides (default 16, recommended 10-16)
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @examples
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#' gr <- GenomicRanges::GRanges(
#'         seqnames = c(PRNP = 'chr20:4699600',             # snp
#'                      HBB  = 'chr11:5227002',             # snp
#'                      HEXA = 'chr15:72346580-72346583',   # del
#'                      CFTR = 'chr7:117559593-117559595'), # ins
#'         strand   = c(PRNP = '+', HBB = '-', HEXA = '-', CFTR = '+'), 
#'         seqinfo  = BSgenome::seqinfo(bsgenome))
#' extend_to_find_gg(gr, plot = TRUE)
#' @export
extend_to_find_gg <- function(gr, nrt=16, plot = FALSE){

    # Extend
    fw <- copy(gr, strand='+', start = end(gr) - nrt + 5, end = start(gr) + 5)
    rv <- copy(gr, strand='-', start = end(gr) - 5, end = start(gr) + nrt - 5)
    fw$targetstart <- rv$targetstart   <- start(gr)
    fw$targetend   <- rv$targetend     <- end(gr)
    fw$targetstrand <- rv$targetstrand <- strand(gr)
    ext <- sort(c(fw, rv))

    # Plot    
    if (plot){
        gr$set <- 'targets'
        fw$set <- 'possible GG of fwd strand'
        rv$set <- 'possible GG on rev strand'
        plot_intervals(c(gr, fw, rv), color_var = 'set', yby = 'set')
    }
    
    # Return
    ext
}

#' Find GG
#' @param gr \code{\link[GenomicRanges]{GRanges-class}}
#' @return   \code{\link[GenomicRanges]{GRanges-class}}
#' @examples
#' require(magrittr)
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#' gr <- GenomicRanges::GRanges(
#'         seqnames = c(PRNP = 'chr20:4699600',             # snp
#'                      HBB  = 'chr11:5227002',             # snp
#'                      HEXA = 'chr15:72346580-72346583',   # del
#'                      CFTR = 'chr7:117559593-117559595'), # ins
#'         strand   = c(PRNP = '+', HBB = '-', HEXA = '-', CFTR = '+'), 
#'         seqinfo  = BSgenome::seqinfo(bsgenome))
#' gr %<>% extend_to_find_gg(plot = TRUE)
#' gr %>% add_seq(bsgenome) %>% find_gg()
#' @export
find_gg <- function(gr){
    
    res <- stri_locate_all_fixed(gr$seq, 'GG')
    gr$substart <- vextract1(res)
    gr$subend   <- vextract2(res)
    
    # Rm GG-free targetranges
    idx <- gr$substart == 'NA'
    if (sum(idx)>0) gr %<>% extract(!idx)

    # Identify GG sites
    gg_dt <- as.data.table(gr) %>% 
            setnames(c('start', 'end'), c('windowstart', 'windowend')) %>% 
            separate_rows(substart, subend)                            %>%
            data.table()                                               %>% 
            extract(, substart := as.numeric(substart))                %>% 
            extract(, subend   := as.numeric(subend))                  %>% 
            # extract(, seq    := substr(seq, substart, subend))  %>%
            extract(strand=='+', start := windowstart + substart - 1)  %>% 
            extract(strand=='+', end   := windowstart + subend   - 1)  %>% 
            extract(strand=='-', end   := windowend   - substart + 1)  %>% 
            extract(strand=='-', start := windowend   - subend   + 1)
    gg_dt %<>% setorderv(c('seqnames', 'start', 'strand'))
    gg_dt[, c('windowstart', 'windowend', 'seq', 'substart', 'subend') := NULL]
    gg_ranges   <- GRanges(gg_dt, seqinfo = seqinfo(gr))
    gg_ranges %<>% add_seq(bsgenome, verbose = FALSE)
    return(gg_ranges)
}


#' Find prime editing sites
#' @param gr        \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome  \code{\link[BSgenome]{BSgenome-class}}
#' @param fix       character vector: fixing sequence on '+' strand
#' @param nprimer   n primer nucleotides (default 13, max 17)
#' @param nrt       n RT nucleotides (default 16, recommended 10-16)
#' @param plot      TRUE (default) or FALSE
#' @examples
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#' gr <- GenomicRanges::GRanges(
#'         seqnames = c(PRNP = 'chr20:4699600',             # snp
#'                      HBB  = 'chr11:5227002',             # snp
#'                      HEXA = 'chr15:72346580-72346583',   # del
#'                      CFTR = 'chr7:117559593-117559595'), # ins
#'         strand   = c(PRNP = '+', HBB = '-', HEXA = '-', CFTR = '+'), 
#'         seqinfo  = BSgenome::seqinfo(bsgenome))
#' BSgenome::getSeq(bsgenome, gr)
#' @export
find_pe_sites <- function(
    gr, 
    bsgenome, 
    fix     = BSgenome::getSeq(bsgenome, gr, as.character = TRUE), 
    nprimer = 13, 
    nrt     = 16, 
    plot    = TRUE
){
    
    # Assert
    assert_is_all_of(gr, 'GRanges')
    assert_is_all_of(bsgenome, 'BSgenome')
    assert_is_character(fix)
    assertive.strings::assert_all_are_matching_regex(fix, '^[ACGTacgt]+$')
    assert_is_a_number(nprimer)
    assert_is_a_number(nrt)
    assertive.numbers::assert_all_are_less_than(nprimer, 17)

    # Find GG in nrt window around target site
    gr$targetname <-  if (has_names(gr)){  names(gr) 
                    } else {               sprintf('target%d', seq_along(gr)) }
    gg <-  extend_to_find_gg(gr, nrt) %>% add_seq(bsgenome) %>% find_gg()
    gg$site <- table(gg$targetname) %>% lapply(seq_len) %>% unlist() %>% names()
    gg$seq <- NULL
    
    # Extract from these other components
    spacerpam <- up_flank(gg, -21,        +1)
    primer    <- up_flank(gg, -5-nprimer, -5)
    fix       <- up_flank(gg, -4,         -4+nrt)
    ext       <- up_flank(gg, -5-nprimer, -4+nrt) %>% invertStrand()
    if (plot){
        spacerpam$part <- 'spacerpam';  primer$part <- 'primer'; fix$part <- 'fix'
        allranges <- c(spacerpam, primer, fix)
        allranges$part %<>% factor(rev(c('spacerpam', 'fix', 'primer')))
        plot_intervals(
            allranges, yby = 'site', color_var = 'part', 
            size_var = 'part', facet_var = c('seqnames', 'targetname'))
        spacerpam$part <- primer$part <- fix$part <- NULL
        
    }
    
    # Return
    spacerpam$spacerpam <- BSgenome::getSeq(bsgenome, spacerpam, as.character = TRUE)
    spacerpam$primer    <- BSgenome::getSeq(bsgenome, primer, as.character = TRUE)
    spacerpam$fix       <- BSgenome::getSeq(bsgenome, fix,    as.character = TRUE)
    spacerpam$extension <- BSgenome::getSeq(bsgenome, ext,    as.character = TRUE)
    spacerpam
}

# PRNP snp: Kuru resistance variant (G -> T)
    # gr  <-  GRanges('chr20:4699500', strand = '+', seqinfo = bsinfo)
    # gr  %<>% multicrispr::add_inverse_strand()
    # (gr %<>% multicrispr::add_seq(bsgenome))
    # (extended <- multicrispr::extend(gr, bsgenome = bsgenome))
    # Biostrings::complement(Biostrings::DNAStringSet(extended$seq[2]))
    # (cas9s <- multicrispr::find_crispr_sites(extended))
    # 
    # # Precision editing
    # #
    # #    ------------^---===....*....
    # # 5' TGGTGGCTGGGG TCAAGGAGGTGGCACCCACAGTCAGTGGAACAA 3'
    # # 3' ACCACCGACCCC AGTTCCTCCACCGTGGGTGTCAGTCACCTTGTT 5'
    # #

# Tay Sachs    
    # gr <-GRanges('chr15:72346580-72346583', strand = '-', 
    #                              seqinfo = seqinfo(bsgenome))
    # 
    # BSgenome::getSeq(bsgenome, gr)
    # gr %<>% multicrispr:::add_inverse_strand()
    # extended <- multicrispr::extend(gr)
    # extended %<>% multicrispr::add_seq(bsgenome)
    # cas9s <- multicrispr::find_crispr_sites(extended)
    
    # Precision editing
    # 
    #                       ===*================---
    # 5' AATGTGAGACAGCTTAAAATAAAATTAACTATAAGAAACTGGTAA
    # 3' TTACCAGTTTCTTATAGTTAATTTTATTTTAAGCTGTCTCACATT
