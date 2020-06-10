
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

#' Extend prime editing target to find GG sites
#' 
#' Extend prime editing target to find GG sites in accessible neighbourhood
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
#' # PE example
#' #-----------
#'     require(magrittr)
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#'     gr <- char_to_granges(c(PRNP = 'chr20:4699600:+',             # snp
#'                            HBB  = 'chr11:5227002:-',              # snp
#'                            HEXA = 'chr15:72346580-72346583:-',    # del
#'                            CFTR = 'chr7:117559593-117559595:+'),  # ins
#'                           bsgenome)
#'     extend_pe_to_gg(gr, plot = TRUE)
#' @export
extend_pe_to_gg <- function(gr, nrt=16, plot = FALSE){

    # Extend
    fw <- copy(gr, strand='+', start = end(gr) - nrt + 5, end = start(gr) + 5)
    rv <- copy(gr, strand='-', start = end(gr) - 5, end = start(gr) + nrt - 5)
    #fw$tstart  <- rv$tstart   <- start(gr)
    #fw$tend    <- rv$tend     <- end(gr)
    #fw$tstrand <- rv$tstrand  <- strand(gr)
    ext <- sort(c(fw, rv))

    # Plot    
    if (plot){
        gr$set <- 'targets'
        fw$set <- 'possible GG of fwd strand'
        rv$set <- 'possible GG on rev strand'
        plot_intervals(c(gr, fw, rv), color_var = 'set', y = 'set')
    }
    
    # Return
    ext
}


#' Find GG
#' @param gr \code{\link[GenomicRanges]{GRanges-class}}
#' @return   \code{\link[GenomicRanges]{GRanges-class}}
#' @examples
#' # PE example
#' #-----------
#'     require(magrittr)
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#'     gr <- char_to_granges(c(PRNP = 'chr20:4699600:+',            # snp
#'                             HBB  = 'chr11:5227002:-',             # snp
#'                             HEXA = 'chr15:72346580-72346583:-',   # del
#'                             CFTR = 'chr7:117559593-117559595:+'), # ins
#'                           bsgenome)
#'     gr %<>% extend_pe_to_gg(plot = TRUE) %>% add_seq(bsgenome) 
#'     find_gg(gr)
#' @export
find_gg <- function(gr){
    
    substart <- subend <- windowstart <- windowend <- NULL
    
    res <- stri_locate_all_fixed(gr$seq, 'GG', overlap=TRUE)
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
            extract(strand=='-', start := windowend   - subend   + 1)  %>% 
            extract()
    gg_dt %<>% setorderv(c('seqnames', 'start', 'strand'))
    gg_dt[, c('windowstart', 'windowend', 'seq', 'substart', 'subend') := NULL]
    gg_ranges   <- GRanges(gg_dt, seqinfo = seqinfo(gr))
    #gg_ranges %<>% add_seq(bsgenome, verbose = FALSE)
    return(gg_ranges)
}


#' Find prime editing spacers
#' 
#' Find prime editing spacers around target ranges
#' 
#' Below the architecture of a prime editing site.
#' Edits can be performed anywhere in the revtranscript area.
#' 
#'         spacer        pam
#'   --------------------===
#'           primer     revtranscript
#'       -------------================
#'   1..............17....GG..........
#'   .....................CC..........
#'       ----------extension----------
#' 
#' @param gr        \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome  \code{\link[BSgenome]{BSgenome-class}}
#' @param edits     character vector: desired edits (on '+' strand)
#' @param nprimer   n primer nucleotides (default 13, max 17)
#' @param nrt       n rev transcr nucleotides (default 16, recomm. 10-16)
#' @param plot      TRUE (default) or FALSE
#' @return  \code{\link[GenomicRanges]{GRanges-class}} with PE spacer ranges
#' Each prime editing range is defined in terms of its N20NGG spacer.
#' Additionally, three sequence mcols are returned:
#'   * spacer: N20 spacers
#'   * pam:    NGG PAMs
#'   * 3pext:  3' extension of gRNA (RTtemplate + primerbindingsite)
#' @examples
#' # Find PE spacers for four clinically relevant loci (Anzalone et al, 2019)
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#'     gr <- char_to_granges(c(
#'         PRNP = 'chr20:4699600:+',             # snp: prion disease
#'         HBB  = 'chr11:5227002:-',             # snp: sickle cell anemia
#'         HEXA = 'chr15:72346580-72346583:-',   # del: tay sachs disease
#'         CFTR = 'chr7:117559593-117559595:+'), # ins: cystic fibrosis
#'         bsgenome)
#'     find_pe_spacers(gr, bsgenome)
#'     find_spacers(extend_for_pe(gr), bsgenome, complement = FALSE)
#'     
#' # Find 
#' @seealso \code{\link{find_spacers}} to find standard crispr sites
#' @export
find_pe_spacers <- function(gr, bsgenome, edits = get_plus_seq(bsgenome, gr), 
    nprimer = 13, nrt = 16, plot = TRUE){
    
    # Assert
    assert_is_all_of(gr, 'GRanges')
    assert_is_all_of(bsgenome, 'BSgenome')
    assert_is_character(edits)
    assert_all_are_matching_regex(edits, '^[ACGTacgt]+$')
    assert_is_a_number(nprimer)
    assert_is_a_number(nrt)
    assert_all_are_less_than(nprimer, 17)

    # Find GG in nrt window around target site
    gr %<>% name_uniquely()
    gg  <-  gr %>% extend_pe_to_gg(nrt) %>% add_seq(bsgenome) %>% find_gg()
    names(gg) <- gg$crisprname <- uniquify(gg$targetname)
    
    # Extract from these other components
    bs <- bsgenome; invstr <- invertStrand
    spacer       <- up_flank(gg, -21,        -2)    %>% add_seq(bs)
    pam          <- up_flank(gg,  -1,        +1)    %>% add_seq(bs) 
    primer       <- up_flank(gg, -4-nprimer, -5)    %>% add_seq(bs)
    revtranscript<- up_flank(gg, -4, -5+nrt) %>% add_fixed_seqs(bs, edits)
    ext  <- up_flank(gg, -4-nprimer, -5+nrt) %>% invertStrand() %>% 
            add_fixed_seqs(bs, edits) # Revcomp for "-" seqs
    # Plot
    if (plot){
        spacer$part<-'spacer'
        ext$part  <- "3' extension"
        allranges <- c(spacer, ext)
        allranges$part %<>% factor((c("spacer", "3' extension")))
        print(plot_intervals(allranges, y = 'crisprname', linetype_var = 'part',
                facet_var = c('seqnames', 'targetname')))
        spacer$part <- NULL
    }
    
    # Add sequences and return
    names(mcols(spacer)) %<>% stri_replace_first_fixed('seq', 'crisprspacer')
    spacer$crisprpam     <- pam$seq
    spacer$primer        <- primer$seq
    spacer$revtranscript <- revtranscript$seq
    spacer$extension     <- ext$seq
    spacer
}

add_fixed_seqs <- function(gr, bsgenome, edits){
    # Get '+' seq
    gr$seq <- get_plus_seq(bsgenome, gr)
    substr(gr$seq, gr$targetstart - start(gr)+1, 
                    gr$targetend  - start(gr)+1) <- edits[gr$targetname]
    gr$seq[as.logical(strand(gr)=='-')] %<>% revcomp()     
    gr
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

