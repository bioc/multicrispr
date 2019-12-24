
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
#' @param plot TRUE or FALSE (default)
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
    
    substart <- subend <- windowstart <- windowend <- NULL
    
    res <- stri_locate_all_fixed(gr$seq, 'GG')
    gr$substart <- vextract1(res)
    gr$subend   <- vextract2(res)
    
    # Rm GG-free gr
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
    #gg_ranges %<>% add_seq(bsgenome, verbose = FALSE)
    return(gg_ranges)
}

name_ranges <- function(gr, name_base = 'range'){
    if (has_names(gr)){  names(gr) 
    } else {             sprintf('%s%d', name_base, seq_along(gr)) }
}


get_plus_seq <- function(bsgenome, gr){
    seqs1  <-  BSgenome::getSeq(bsgenome, 
                                seqnames(gr),
                                start(gr),
                                end(gr), 
                                strand = '+', 
                                as.character = TRUE)
    if (has_names(gr)) names(seqs1) <- names(gr)
    seqs1
        
}


revcomp <- function(y)  y %>% 
                        Biostrings::DNAStringSet() %>% 
                        Biostrings::reverseComplement() %>% 
                        as.character()


#' Find prime editing sites
#' 
#' Find prime editing sites around target ranges
#' 
#' Below the architecture of a prime editing site.
#' Fixes can be performed anywhere in the revtranscription area.
#' 
#'           primer   revtranscription
#'       -------------================
#'   1..............17....GG..........
#'   --------------------===
#'         spacer        pam
#' 
#' @param gr        \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome  \code{\link[BSgenome]{BSgenome-class}}
#' @param fixes     character vector: '+' strand fix seqs
#' @param nprimer   n primer nucleotides (default 13, max 17)
#' @param nrt       n rev transcr nucleotides (default 16, recomm. 10-16)
#' @param plot      TRUE (default) or FALSE
#' @return  \code{\link[GenomicRanges]{GRanges-class}} with prime editing sites.
#' Each prime editing range is defined in terms of its N20NGG spacer.
#' Additionally, five sequences are returned per PE site:
#'   * spacer: N20 spacers
#'   * pam:    NGG PAMs
#'   * primer: primers (end at N17 of spacer)
#'   * revtrans: reverse transcription sequences (start at N18 of spacer)
#'   * ext:    3' extensions of gRNA (reverse complements of primer+revtrans)
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
#' find_pe_sites(gr, bsgenome)
#' @export
find_pe_sites <- function(gr, bsgenome, fixes = get_plus_seq(bsgenome, gr), 
    nprimer = 13, nrt = 16, plot = TRUE){
    
    # Assert
    assert_is_all_of(gr, 'GRanges')
    assert_is_all_of(bsgenome, 'BSgenome')
    assert_is_character(fixes)
    assert_all_are_matching_regex(fixes, '^[ACGTacgt]+$')
    assert_is_a_number(nprimer)
    assert_is_a_number(nrt)
    assert_all_are_less_than(nprimer, 17)

    # Find GG in nrt window around target site
    gr$targetname <- name_ranges(gr, 'target')
    gg  <-  gr %>% extend_to_find_gg(nrt) %>% add_seq(bsgenome) %>% find_gg()
    gg$site <- table(gg$targetname) %>% lapply(seq_len) %>% unlist() %>% names()
    gg$seq <- NULL
    
    # Extract from these other components
    spacer   <- up_flank(gg, -21,        -2)     %>% add_seq(bsgenome)
    pam      <- up_flank(gg,  -1,        +1)     %>% add_seq(bsgenome) 
    primer   <- up_flank(gg, -4-nprimer, -5)     %>% add_seq(bsgenome)
    revtrans <- up_flank(gg, -4,         -5+nrt) %>% add_seq(bsgenome)
    ext      <- up_flank(gg, -4-nprimer, -5+nrt) %>% invertStrand()
    ext$seq  <- get_plus_seq(bsgenome, ext)              # Get "+" seq
    substr(ext$seq, ext$targetstart-start(ext)+1, 
                    ext$targetend  -start(ext)+1) <- fixes[ext$targetname]
    ext$seq[as.logical(strand(ext)=='-')] %<>% revcomp() # Revcomp for "-" seqs

    # Plot
    if (plot){
        spacer$part<-'spacer'; primer$part<-'primer'; revtrans$part<-'revtrans'
        allranges <- c(spacer, primer, revtrans)
        allranges$part %<>% factor(rev(c('spacer', 'revtrans', 'primer')))
        plot_intervals(allranges, yby = 'site', color_var = 'part', 
            size_var = 'part', facet_var = c('seqnames', 'targetname'))
        spacer$part <- primer$part <- revtrans$part <- NULL
    }
    
    # Add sequences and return
    names(mcols(spacer)) %<>% stri_replace_first_fixed('seq', 'spacer')
    spacer$pam <- pam$seq 
    spacer$primer<-primer$seq; 
    spacer$revtrans<-revtrans$seq
    spacer$ext <- ext$seq
    spacer
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
