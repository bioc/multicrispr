
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



add_nickspacers <- function(primespacers, bsgenome, 
    ontargets = c('Doench2014', 'Doench2016'), mismatches = 2, 
    offtargetmethod = c('bowtie', 'vcountpdict')[1],
    indexedgenomesdir = INDEXEDGENOMESDIR, outdir = OUTDIR, plot = TRUE
){
# Check / Initialize
    . <- crisprname <- crisprspacer <- crisprpam <- NULL
    off0 <- off1 <- off2 <- pename <- off <- NULL
    primespacers$pename <- primespacers$crisprname
    primespacers$type <- 'pespacer'
# Get nick spacers
    message('\tAdd nickspacers')
    nickzone <- invertStrand(down_flank(primespacers, -3+40-5, -3+90+17))
    mcols(nickzone) %<>% extract(, c('targetname', 'pename'), drop = FALSE)
    nickspacers <- find_spacers(nickzone, bsgenome, complement = FALSE, 
        ontargets = ontargets, offtargetmethod = offtargetmethod, 
        offtargetfilterby = 'pename',
        mismatches = mismatches, indexedgenomesdir = indexedgenomesdir, 
        outdir = outdir, plot=FALSE)
    nickspacers$type <- 'nickspacer'
    if (verbose) message('\tFound ', length(nickspacers), ' nickspacers')
    mcols(nickspacers)[[ontargets]] %<>% round(digits = 2)
# Merge
    nickdt  <-  gr2dt(nickspacers) %>% 
                extract( , .(
                    pename     = pename,
                    nickrange  = as.character(granges(nickspacers)), 
                    nickspacer = crisprspacer, 
                    nickpam    = crisprpam,
                    nickoff    = off,
                    nickoff0   = off0, 
                    nickoff1   = off1, 
                    nickoff2   = off2)) %>%
                extract(, (paste0('nick', ontargets)) := 
                                mcols(nickspacers)[[ontargets]]) %>% 
                extract( , lapply(.SD, pastelapse), by = 'pename')

    pedt <- gr2dt(primespacers)
    pedt$type <- NULL
    mergedranges <- merge(pedt, nickdt, by = 'pename', sort = FALSE, all = TRUE)
    mergedranges$pename <- NULL
    mergedranges %<>% dt2gr(seqinfo(primespacers))
# Plot and Return
    if (plot)   print(plot_intervals(mergedranges))
    return(mergedranges)
}


pastelapse <- function(x) paste0(x, collapse = ';')

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
#' @param edits     character vector: desired edits on '+' strand.
#'                  If named, names should be identical to those of \code{gr}
#' @param nprimer   n primer nucleotides (default 13, max 17)
#' @param nrt       n rev transcr nucleotides (default 16, recomm. 10-16)
#' @param ontargets  'Doench2014' or 'Doench2016': on-target scoring method
#' @param mismatches  no of primespacer mismatches 
#'                   (default 0, to suppress offtarget analysis: -1)
#' @param nickmatches no of nickspacer offtarget mismatches 
#'                   (default 2, to suppresses offtarget analysis: -1)
#' @param offtargetmethod  'bowtie' or 'vcountpdict'
#' @param indexedgenomesdir  directory with indexed genomes 
#'                           (as created by \code{\link{index_genome}})
#' @param outdir    directory whre offtarget analysis output is written
#' @param verbose   TRUE (default) or FALSE
#' @param plot      TRUE (default) or FALSE
#' @param ...       passed to plot_intervals
#' @return  \code{\link[GenomicRanges]{GRanges-class}} with prime editing spacer
#' ranges and following mcols: 
#'   * crisprspacer: N20 spacers
#'   * crisprpam:    NGG PAMs
#'   * crisprprimer: primer (on PAM strand)
#'   * crisprtranscript: reverse transcript (on PAM strand)
#'   * crisprextension:  3' extension of gRNA
#'               contains: reverse transcription template + primer binding site
#'               sequence can be found on non-PAM strand
#'   * crisprextrange: genomic range of crispr extension
#'   * Doench2016|4: on-target efficiency scores
#'   * off0, off1, off2: number of offtargets with 0, 1, 2 mismatches
#'   * off: total number of offtargets: off = off0 + off1 + ...
#'   * nickrange:        nickspacer range
#'   * nickspacer:       nickspacer sequence
#'   * nickDoench2016|4: nickspacer Doench scores
#'   * nickoff:          nickspacer offtarget counts
#' @examples
#' # Find PE spacers for 4 clinically relevant loci (Anzalone et al, 2019)
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#'     gr <- char_to_granges(c(
#'         PRNP = 'chr20:4699600:+',             # snp: prion disease
#'         HBB  = 'chr11:5227002:-',             # snp: sickle cell anemia
#'         HEXA = 'chr15:72346580-72346583:-',   # del: tay sachs disease
#'         CFTR = 'chr7:117559593-117559595:+'), # ins: cystic fibrosis
#'         bsgenome)
#'     spacers <- find_primespacers(gr, bsgenome)
#'     spacers <- find_spacers(extend_for_pe(gr), bsgenome, complement = FALSE)
#'     
#' # Edit PRNP locus for resistance against prion disease (Anzalone et al, 2019)
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#'     gr <- char_to_granges(c(PRNP = 'chr20:4699600:+'), bsgenome)
#'     find_primespacers(gr, bsgenome)
#'     find_primespacers(gr, bsgenome, edits = 'T')
#' @seealso \code{\link{find_spacers}} to find standard crispr sites
#' @export
find_primespacers <- function(gr, bsgenome, edits = get_plus_seq(bsgenome, gr), 
    nprimer = 13, nrt = 16, ontargets = c('Doench2014', 'Doench2016')[1], 
    mismatches = 0, nickmatches = 2, 
    offtargetmethod = c('bowtie', 'vcountpdict')[1], 
    indexedgenomesdir = INDEXEDGENOMESDIR, outdir = OUTDIR, 
    verbose = TRUE, plot = TRUE, ...){
# Assert
    assert_is_all_of(gr, 'GRanges')
    assert_is_all_of(bsgenome, 'BSgenome')
    assert_is_character(edits)
    assert_all_are_matching_regex(edits, '^[ACGTacgt]+$')
    assert_are_same_length(gr, edits)
    if (has_names(edits))   assert_are_identical(names(gr), names(edits))
    assert_is_a_number(nprimer)
    assert_is_a_number(nrt)
    assert_all_are_less_than(nprimer, 17)
    assert_is_subset(ontargets, c('Doench2016', 'Doench2014'))
    assert_is_a_bool(verbose)
    assert_is_a_bool(plot)
# Find GG in nrt window around target site
    gr %<>% name_uniquely(); names(edits) <- names(gr)
    gg  <-  gr %>% extend_pe_to_gg(nrt) %>% add_seq(bsgenome) %>% find_gg()
    names(gg) <- gg$crisprname <- uniquify(gg$targetname)
# Extract from these other components
    bs <- bsgenome; invstr <- invertStrand
    spacers       <- up_flank(gg, -21,        -2)    %>% add_seq(bs)
    pams          <- up_flank(gg,  -1,        +1)    %>% add_seq(bs) 
    primers       <- up_flank(gg, -4-nprimer, -5)    %>% add_seq(bs)
    transcripts   <- up_flank(gg, -4, -5+nrt) %>% add_fixed_seqs(bs, edits)
    exts  <- up_flank(gg, -4-nprimer, -5+nrt) %>% invertStrand() %>% 
            add_fixed_seqs(bs, edits) # Revcomp for "-" seqs
# Add sequences
    names(mcols(spacers)) %<>% stri_replace_first_fixed('seq', 'crisprspacer')
    spacers$crisprpam        <- pams$seq
    spacers$crisprprimer     <- primers$seq
    spacers$crisprtranscript <- transcripts$seq
    spacers$crisprextension  <- exts$seq
    spacers$crisprextrange   <- unname(as.character(granges(exts)))
# Add offtargets, ontargets, nickspacers
    if (verbose) message("\tFound ", length(spacers), " spacers/3'extensions")
    spacers %<>% add_offtargets(bsgenome, mismatches = mismatches, pam = 'NGG',
                    offtargetmethod = offtargetmethod,outdir = outdir, 
                    indexedgenomesdir = indexedgenomesdir, 
                    verbose= TRUE, plot = FALSE)
    spacers %<>% add_ontargets(bsgenome, method = ontargets, plot = FALSE)
    spacers %<>% add_nickspacers(bsgenome, ontargets = ontargets, 
            mismatches = nickmatches, offtargetmethod = offtargetmethod,
            indexedgenomesdir = indexedgenomesdir,
            outdir = outdir, plot = FALSE)
# Plot and Return
    if (plot) print(plot_intervals(spacers, ...))
    spacers
}

find_pe_spacers <- function(...){
    .Deprecated('find_primespacers')
    find_primespacers(...)
}

add_fixed_seqs <- function(gr, bsgenome, edits){
    # Get '+' seq
    gr$seq <- get_plus_seq(bsgenome, gr)
    substr(gr$seq, gr$targetstart - start(gr)+1, 
                    gr$targetend  - start(gr)+1) <- edits[gr$targetname]
    gr$seq[as.logical(strand(gr)=='-')] %<>% revcomp()     
    gr
}

