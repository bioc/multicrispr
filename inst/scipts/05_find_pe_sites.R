
cextract1 <- function(y) y[, 1] %>% paste0(collapse=';')
cextract2 <- function(y) y[, 2] %>% paste0(collapse=';')
vextract1 <- function(z) vapply(z, cextract1, character(1))
vextract2 <- function(z) vapply(z, cextract2, character(1))

#' Find prime editing sites
#' @param gr        \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome  \code{\link[BSgenome]{BSgenome-class}}
#' @param fix       character vector: fixing sequence on '+' strand
#' @param npb       number (default 13, max 17) of primer binding nucleotides
#' @param nrt       number (default 16) of reverse transcription nucleotides. Anzalone et al. (2019) recommend 10-16.
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
#' fix <- c('T', 'A')
#' @export
find_pe_sites <- function(gr, bsgenome, fix, npb = 13, nrt = 16){
    
    # Assert
    assert_is_all_of(gr, 'GRanges')
    assert_is_all_of(bsgenome, 'BSgenome')
    assert_is_character(fix)
    assertive.strings::assert_all_are_matching_regex(fix, '^[ACGTacgt]+$')
    assert_is_a_number(npb)
    assert_is_a_number(nrt)
    assertive.numbers::assert_all_are_less_than(npb, 17)

    # Find GG in revtrans window around target site
    gr$targetname <- if (has_names(gr))  names(gr) else  sprintf('TRG%d', seq_along(gr))
    gr %<>% add_seq(bsgenome)
    fw <- rv <- gr
    fw$targetstart  <- rv$targetstart  <- start(gr)
    fw$targetend    <- rv$targetend    <- end(gr)
    fw$targetstrand <- rv$targetstrand <- strand(gr)
    strand(fw) <- '+'; start(fw) <- end(gr)-nrt+5;  end(fw) <- start(gr)+5
    strand(rv) <- '-'; start(rv) <- end(gr)-5;      end(rv) <- start(gr)+nrt-5
    ggwindow <- c(fw, rv) %>% sort() %>% add_seq(bsgenome)
    res <- stri_locate_all_fixed(ggwindow$seq, 'GG')
    ggwindow$substart <- vextract1(res)
    ggwindow$subend   <- vextract2(res)
    
    # Rm GG-free targetranges
    idx <- ggwindow$substart == 'NA'
    if (sum(idx)>0) ggwindow %<>% extract(!idx)

    # Identify GG sites
    sites_dt <- as.data.table(ggwindow) %>% 
                data.table::setnames(c('start', 'end'), c('windowstart', 'windowend')) %>% 
                tidyr::separate_rows(substart, subend)                %>%
                data.table()                                          %>% 
                extract(, substart := as.numeric(substart))           %>% 
                extract(, subend   := as.numeric(subend))             %>% 
                # extract(, seq    := substr(seq, substart, subend))  %>%
                extract(strand=='+', start := windowstart + substart - 1) %>% 
                extract(strand=='+', end   := windowstart + subend   - 1) %>% 
                extract(strand=='-', end   := windowend   - substart + 1) %>% 
                extract(strand=='-', start := windowend   - subend   + 1)
    sites_dt %<>% setorderv(c('seqnames', 'start', 'strand'))
    sites_dt %>%  extract(, ggname := sprintf('%s_%d', targetname, 1:.N), 
                            by = 'targetname')
    sites_dt[, c('windowstart', 'windowend', 'seq', 'substart', 'subend') := NULL]
    ggsites   <- GRanges(sites_dt, seqinfo = seqinfo(gr))
    ggsites %<>% add_seq(bsgenome)
    inv <- invertStrand
    spacer   <-ggsites %>% up_flank(-21,     -2,     y_by = 'ggname', bsgenome = bsgenome)
    pbs      <- spacer %>% up_flank( 17-npb, 17,     y_by = 'ggname', bsgenome = bsgenome) %>% inv()
    template <- spacer %>% down_flank( -3,     -3+nrt, y_by = 'ggname', bsgenome = bsgenome) %>% inv()
    extension<- spacer %>%   extend( 17-npb, -3+nrt, y_by = 'ggname', bsgenome = bsgenome) %>% inv()
    
    spacer$component    <- 'spacer'
    pbs$component       <- 'primerbinding'
    template$component  <- 'revtranscription'
    gr$component   <- 'target'
    allranges <- c(spacer, pbs, template)
    allranges$component %<>% factor(rev(c('spacer', 'revtranscription', 'primerbinding')))
    plot_intervals(allranges, y_by = 'ggname', color_var = 'component', size_var = 'component', facet_var = c('seqnames', 'targetname'))


}

# PRNP snp: Kuru resistance variant (G -> T)
    gr  <-  GRanges('chr20:4699500', strand = '+', seqinfo = bsinfo)
    gr  %<>% multicrispr::add_inverse_strand()
    (gr %<>% multicrispr::add_seq(bsgenome))
    (extended <- multicrispr::extend(gr, bsgenome = bsgenome))
    Biostrings::complement(Biostrings::DNAStringSet(extended$seq[2]))
    (cas9s <- multicrispr::find_crispr_sites(extended))

    # Precision editing
    #
    #    ------------^---===....*....
    # 5' TGGTGGCTGGGG TCAAGGAGGTGGCACCCACAGTCAGTGGAACAA 3'
    # 3' ACCACCGACCCC AGTTCCTCCACCGTGGGTGTCAGTCACCTTGTT 5'
    #

# Tay Sachs    
    gr <-GRanges('chr15:72346580-72346583', strand = '-', 
                                 seqinfo = seqinfo(bsgenome))
    
    BSgenome::getSeq(bsgenome, gr)
    gr %<>% multicrispr:::add_inverse_strand()
    extended <- multicrispr::extend(gr)
    extended %<>% multicrispr::add_seq(bsgenome)
    cas9s <- multicrispr::find_crispr_sites(extended)
    
    # Precision editing
    # 
    #                       ===*================---
    # 5' AATGTGAGACAGCTTAAAATAAAATTAACTATAAGAAACTGGTAA
    # 3' TTACCAGTTTCTTATAGTTAATTTTATTTTAAGCTGTCTCACATT
