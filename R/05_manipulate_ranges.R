#' Add sequence to GRanges
#' @param gr           \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome     \code{\link[BSgenome]{BSgenome-class}}
#' @param verbose      TRUE or FALSE (default)
#' @param as.character TRUE (default) or FALSE
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @examples 
#' # PE example
#' #-----------
#'     require(magrittr)
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#'     gr <- char_to_granges(c(PRNP = 'chr20:4699600:+',             # snp
#'                             HBB  = 'chr11:5227002:-',             # snp
#'                             HEXA = 'chr15:72346580-72346583:-',   # del
#'                             CFTR = 'chr7:117559593-117559595:+'), # ins
#'                           bsgenome)
#'    (gr %<>% add_seq(bsgenome))
#'    
#' # TFBS example
#' #-------------
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'     bedfile  <- system.file('extdata/SRF.bed', package='multicrispr')
#'     gr <- bed_to_granges(bedfile, 'mm10')
#'     (gr %<>% add_seq(bsgenome))
#' @export
add_seq <- function(gr, bsgenome, verbose = FALSE, as.character = TRUE){
    
    # Assert
    assert_is_all_of(gr, 'GRanges')
    assert_is_all_of(bsgenome, 'BSgenome')
    
    # Message
    if (verbose)  cmessage('\tAdd seq')
    
    # Align seqlevelsStyle if required
    if (seqlevelsStyle(bsgenome)[1] != seqlevelsStyle(gr)[1]){
            cmessage("\t\t\tSet seqlevelStyle(bsgenome) <- seqlevelStyle(gr)")
            seqlevelsStyle(bsgenome)[1] <- seqlevelsStyle(gr)[1]
    }
    
    # Add seq
    gr$seq <- unname(BSgenome::getSeq(
                        bsgenome,
                        names        = seqnames(gr),
                        start        = start(gr),
                        end          = end(gr), 
                        strand       = strand(gr), 
                        as.character = as.character))
    
    # Return
    gr
}


summarize_loci <- function(gr){
    sprintf('%s:%s-%s', 
            as.character(seqnames(gr)), 
            start(gr), 
            end(gr))
}

#' Extend or Flank GRanges
#' 
#' Returns extensions, upstream flanks, or downstream flanks
#' 
#' \code{up_flank}   returns upstream flanks, in relation to start(gr).
#' \code{down_flank} returns downstream flanks, in relation to end(gr).
#' \code{extend}     returns extensions, in relation to start(gr) and end(gr)
#' @param gr        \code{\link[GenomicRanges]{GRanges-class}}
#' @param start     (pos or neg) number: relative start position (see details)
#' @param end       (pos or neg) number: relative end position   (see details)
#' @param strandaware  TRUE (default) or FALSE: consider strand information?
#' @param bsgenome  NULL (default) or \code{\link[BSgenome]{BSgenome-class}}.
#'                  Required to update gr$seq if present.
#' @param verbose   TRUE or FALSE (default)
#' @param plot      TRUE or FALSE (default)
#' @param linetype_var string: gr var mapped to linetype
#' @param ...       passed to \code{\link{plot_intervals}}
#' @return a \code{\link[GenomicRanges]{GRanges-class}}
#' @examples 
#' # PE example
#' #-----------
#' require(magrittr)
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#' gr <- char_to_granges(c(PRNP  = 'chr20:4699600:+',         # snp
#'                          HBB  = 'chr11:5227002:-',            # snp
#'                          HEXA = 'chr15:72346580-72346583:-',  # del
#'                          CFTR = 'chr7:117559593-117559595:+'),# ins
#'                       bsgenome = bsgenome)
#' gr %>% up_flank( -22,  -1, plot=TRUE, facet_var=c('targetname', 'seqnames'))
#' gr %>% up_flank( -22,  -1, plot=TRUE, strandaware=FALSE)
#' gr %>% down_flank(+1, +22, plot=TRUE)
#' gr %>% down_flank(+1, +22, plot=TRUE, strandaware=FALSE)
#' gr %>% extend(   -10, +20, plot=TRUE)
#' gr %>% extend(   -10, +20, plot=TRUE, strandaware=FALSE)
#'
#' # TFBS example
#' #-------------
#'     bedfile <- system.file('extdata/SRF.bed', package='multicrispr')
#'     gr <- bed_to_granges(bedfile, genome = 'mm10')
#'     gr %>% extend(plot = TRUE)
#'     gr %>% up_flank(plot = TRUE)
#'     gr %>% down_flank(plot = TRUE)
#' @export
up_flank <- function(
    gr, start = -200, end = -1, strandaware = TRUE, bsgenome = NULL, 
    verbose = FALSE, plot = FALSE, linetype_var = 'set', ...
){
    # Assert
    assert_is_all_of(gr, 'GRanges')
    assert_is_a_number(start)
    assert_is_a_number(end)
    assert_is_a_bool(verbose)
    
    # Record
    newgr <- gr
    shift <- sprintf('(%s%d,%s%d)', 
                                csign(start), abs(start), csign(end), abs(end))
    txt <- sprintf('\t\t%d%supstream %s flanks', 
                    length(newgr), 
                    ifelse(!strandaware, ' (strandagnostic) ', ' '),
                    shift)
    
    # Flank
    GenomicRanges::start(newgr) <- GenomicRanges::start(gr) + start
    GenomicRanges::end(newgr)   <- GenomicRanges::start(gr) + end
    if (strandaware){
        idx <- as.logical(strand(newgr)=='-')
        # do not switch following lines to avoid integrity errors !
        GenomicRanges::end(  newgr)[idx] <- GenomicRanges::end(gr)[idx] - start
        GenomicRanges::start(newgr)[idx] <- GenomicRanges::end(gr)[idx] - end
    }

    # Add seq
    if ('seq' %in% names(mcols(gr))){
        assert_is_all_of(bsgenome, 'BSgenome')
        newgr %<>% add_seq(bsgenome)
    }

    # Plot, Message, Return
    if (plot){
        gr$set    <- 'original'
        newgr$set <- 'upstream flanks'
        allgr <- c(gr, newgr)
        allgr$set %<>% factor(c('original', 'upstream flanks'))
        print(plot_intervals(
                allgr, linetype_var = linetype_var, ..., title = txt))
        newgr$set <- NULL
    }
    if (verbose) message(txt)
    newgr
}


#' @rdname up_flank
#' @export
down_flank <- function(
    gr, start = 1,  end = 200, strandaware = TRUE, bsgenome = NULL,
    verbose = FALSE, plot = FALSE, linetype_var = 'set', ...
){
    # Assert
    assert_is_any_of(gr, 'GRanges')
    assert_is_a_number(start)
    assert_is_a_number(end)
    assert_is_a_bool(verbose)
    
    # Record 
    newgr <- gr
    shift <- sprintf('(%s%d,%s%d)', 
                                csign(start), abs(start), csign(end), abs(end))
    txt <- sprintf('\t\t%d%sdownstream %s flanks', 
                    length(newgr), 
                    ifelse(!strandaware, ' (strandagnostic) ', ' '),
                    shift)
    
    # Flank
    GenomicRanges::end(newgr)   <- GenomicRanges::end(gr) + end
    GenomicRanges::start(newgr) <- GenomicRanges::end(gr) + start
    if (strandaware){
        idx <- as.logical(strand(newgr)=='-')
        GenomicRanges::start(newgr)[idx] <- GenomicRanges::start(gr)[idx]-end
        GenomicRanges::end(  newgr)[idx] <- GenomicRanges::start(gr)[idx]-start
    }

    # Add seq
    if ('seq' %in% names(mcols(gr))){
        assert_is_all_of(bsgenome, 'BSgenome')
        newgr %<>% add_seq(bsgenome)
    }
    
    # Plot, Message, Return
    if (plot){
        gr$set <- 'original'
        newgr$set <- 'downstream flanks'
        allgr <- c(gr, newgr)
        allgr$set %<>% factor(c('original', 'downstream flanks'))
        print(plot_intervals(
                allgr, linetype_var = linetype_var, ..., title=txt))
        newgr$set <- NULL
    }
    if (verbose) message(txt)
    newgr
}


#' @rdname up_flank
#' @export
extend <- function(
    gr, start = -22, end = 22, strandaware = TRUE, bsgenome = NULL,
    verbose = FALSE, plot = FALSE, linetype_var = 'set', ...
){

    # Assert
    assert_is_any_of(gr, 'GRanges')
    assert_is_a_number(start)
    assert_is_a_number(end)
    assert_is_a_bool(verbose)
    
    # Record
    newgr <- gr
    shift <- sprintf('(%s%d,%s%d)', 
                    csign(start), abs(start), csign(end), abs(end))
    txt <- sprintf('\t\t%d%sextensions %s', 
                    length(newgr), 
                    ifelse(!strandaware, ' (strandagnostic) ', ' '),
                    shift)

    # Extend
    GenomicRanges::start(newgr) <- GenomicRanges::start(newgr) + start
    GenomicRanges::end(  newgr) <- GenomicRanges::end(  newgr) + end
    if (strandaware){
        idx <- as.logical(strand(newgr)=='-')
        GenomicRanges::start(newgr)[idx] <- GenomicRanges::start(gr)[idx] - end
        GenomicRanges::end(  newgr)[idx] <- GenomicRanges::end(gr)[idx] - start
    }

    # Add seq
    if ('seq' %in% names(mcols(gr))){
        assert_is_all_of(bsgenome, 'BSgenome')
        newgr %<>% add_seq(bsgenome)
    }
    
    # Plot, Message, Return
    if (plot){
        gr$set <- 'original'
        newgr$set <- 'extensions'
        allgr <- c(gr, newgr)
        allgr$set %<>% factor(c('original', 'extensions'))
        print(plot_intervals(
                allgr, linetype_var = linetype_var, ..., title=txt))
        newgr$set <- NULL
    }
    if (verbose) message(txt)
    newgr
    
}




#' Add inverse strand
#' @param gr         \code{\link[GenomicRanges]{GRanges-class}}
#' @param verbose    TRUE or FALSE (default)
#' @param plot       TRUE or FALSE (default)
#' @param ...         \code{\link{plot_intervals}} arguments
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @examples
#' # PE example
#' #-----------
#'     require(magrittr)
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#'     gr <- char_to_granges(c(PRNP = 'chr20:4699600:+',             # snp
#'                             HBB  = 'chr11:5227002:-',             # snp
#'                             HEXA = 'chr15:72346580-72346583:-',   # del
#'                             CFTR = 'chr7:117559593-117559595:+'), # ins
#'                           bsgenome)
#'     add_inverse_strand(gr, plot = TRUE)
#' # TFBS example
#' #-------------
#'     bedfile <- system.file('extdata/SRF.bed', package='multicrispr')
#'     gr <- bed_to_granges(bedfile, genome = 'mm10')
#'     add_inverse_strand(gr)
#' @export
add_inverse_strand <- function(gr, verbose = FALSE, plot = FALSE, ...){

    # Assert
    assertive::assert_is_all_of(gr, 'GRanges')

    # Invert
    revcomps <- invertStrand(gr)
    if ('seq' %in% names(mcols(gr))){
        revcomps$seq <- as.character(Biostrings::reverseComplement(
                            DNAStringSet(gr$seq)))
    }
    
    # Concatenate
    newgr <- c(gr, revcomps)
    newgr %<>% unique()
    txt <- sprintf('\t\t%d ranges after adding inverse strands', length(newgr))
    
    # Sort
    newgr <- sortSeqlevels(newgr)
    newgr <- GenomicRanges::sort(newgr, ignore.strand = TRUE)
    names(newgr) %<>% paste0('_', 
                            c('+'='f', '-'='r')[as.character(strand(newgr))])
    
    # Plot
    if (plot){
        gr$set    <- 'original'
        revcomps$set <- 'inverse'
        print(plot_intervals(
                c(gr, revcomps), 
                color_var = 'set', linetype_var = 'set', ..., title = txt))
        gr$set <- NULL
    }
    
    # Message
    if (verbose) cmessage(txt)
    newgr
}

#' Double flank
#' 
#' @param gr \code{\link[GenomicRanges]{GRanges-class}}
#' @param upstart      upstream flank start in relation to start(gr)
#' @param upend        upstream flank end   in relation to start(gr)
#' @param downstart    downstream flank start in relation to end(gr)
#' @param downend      downstream flank end   in relation to end(gr)
#' @param strandaware  TRUE (default) or FALSE
#' @param plot         TRUE or FALSE (default)
#' @param linetype_var gr var mapped to linetype
#' @param ...          passed to plot_intervals
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @examples 
#' # Prime Editing example
#' #----------------------
#'     require(magrittr)
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#'     gr <- char_to_granges(c(PRNP = 'chr20:4699600:+',             # snp
#'                             HBB  = 'chr11:5227002:-',             # snp
#'                             HEXA = 'chr15:72346580-72346583:-',   # del
#'                             CFTR = 'chr7:117559593-117559595:+'), # ins
#'                           bsgenome)
#'     double_flank(gr, -10,  -1, +1, +20, plot = TRUE)
#'       
#' # TFBS example
#' #-------------
#'     bedfile  <- system.file('extdata/SRF.bed', package='multicrispr')
#'     gr  <-  bed_to_granges(bedfile, genome = 'mm10', plot = FALSE)
#'     double_flank(gr, plot = TRUE)
#' @export
double_flank <- function(
    gr, 
    upstart     = -200, 
    upend       = -1, 
    downstart   = 1, 
    downend     = 200,
    strandaware = TRUE,
    plot        = FALSE, 
    linetype_var = 'set',
    ...
){

    # Up flank, down flank, concatenate
    up <- up_flank(gr,   upstart,   upend,   strandaware, verbose = FALSE)
    dn <- down_flank(gr, downstart, downend, strandaware, verbose = FALSE)
    names(up) %<>% paste0('_u') # ensure unique names
    names(dn) %<>% paste0('_d')
    newgr  <- c(up, dn)
    txt <- sprintf('\t\t%d flank ranges: %d up + %d down', 
                        length(newgr), length(up), length(dn))
    
    # Plot    
    if (plot){
        gr$set <- 'original'
        up$set <- 'upstream flank'
        dn$set <- 'downstream flank'
        allgr <- c(gr, up, dn)
        allgr$set %<>% factor(
                        c('original', 'upstream flank', 'downstream flank'))
        print(plot_intervals(allgr, linetype_var = linetype_var, title = txt, 
                            y = 'targetname', ...))
    }
    
    # Return
    message(txt)
    newgr
}


