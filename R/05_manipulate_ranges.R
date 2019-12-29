#' Add sequence to GRanges
#' @param gr           \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome     \code{\link[BSgenome]{BSgenome-class}}
#' @param verbose      TRUE or FALSE (default)
#' @param as.character TRUE (default) or FALSE
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @examples 
#' # SRF binding sites
#'     bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'     gr <- bed_to_granges(bedfile, 'mm10')
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'     (gr <- add_seq(gr, bsgenome))
#'     
#' # PRNP
#'     require(magrittr)
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#'     gr  <-  GenomicRanges::GRanges(
#'                 'chr20:4699500', strand = '+', 
#'                  seqinfo = BSgenome::seqinfo(bsgenome))
#'     gr %<>% multicrispr::add_inverse_strand()
#'     gr %<>% multicrispr::extend(bsgenome = bsgenome)
#'    (gr %<>% multicrispr::add_seq(bsgenome))
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
#' Operations are performed in a strand-aware fashion by default, but can be 
#' made strand-agnostic by setting stranded=FALSE.
#' 
#' Return values (note: if stranded=FALSE, only + versions are returned)
#'     up_flank  grstart + start  -->  grstart + end   (+)
#'               grend - end     <--   grend - start   (-)
#'            
#'   down_flank  grend + start    -->  grend + end     (+) 
#'               grstart - end    <--  grstart - start (-)
#'               
#'       extend  grstart + start  -->  grend + end     (+)
#'               grstart - end    <--  grend - start   (-)
#' @param gr        \code{\link[GenomicRanges]{GRanges-class}}
#' @param start     (pos or neg) number: relative start position (see details)
#' @param end       (pos or neg) number: relative end position   (see details)
#' @param stranded  TRUE (default) or FALSE: consider strand information
#' @param bsgenome  NULL (default) or \code{\link[BSgenome]{BSgenome-class}}.
#'                  Required to update gr$seq if present.
#' @param verbose   TRUE or FALSE (default)
#' @param plot      TRUE or FALSE (default)
#' @param ...       passed to \code{\link{plot_intervals}}
#' @return a \code{\link[GenomicRanges]{GRanges-class}}
#' @examples 
#' require(magrittr)
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#' gr <- GenomicRanges::GRanges(
#'           c(locus1 = 'chr11:5227002', locus2 = 'chr11:4699500'), 
#'             strand = c('-', '+'), 
#'             seqinfo = BSgenome::seqinfo(bsgenome))
#' gr %>%   up_flank(-22,  -1, plot = TRUE)
#' gr %>%   up_flank(-22,  -1, plot = TRUE, stranded = FALSE)
#' gr %>% down_flank( +1, +22, plot = TRUE)
#' gr %>% down_flank( +1, +22, plot = TRUE, stranded = FALSE)
#' gr %>%     extend(-10, +20, plot = TRUE)
#' gr %>%     extend(-10, +20, plot = TRUE, stranded = FALSE)
#' @export
up_flank <- function(
  gr, 
  start    = -200,
  end      = -1,
  stranded = TRUE,
  bsgenome = NULL,
  verbose  = FALSE,
  plot     = FALSE,
  ...
){
    # Assert
    assert_is_any_of(gr, 'GRanges')
    assert_is_a_number(start)
    assert_is_a_number(end)
    assert_is_a_bool(verbose)
    
    # Left flank
    newgr <- gr
    GenomicRanges::start(newgr) <- GenomicRanges::start(gr) + start
    GenomicRanges::end(newgr)   <- GenomicRanges::start(gr) + end
    formatstr <- '\t\t%d left  flanks: [start%s%d, start%s%d]'
    if (stranded){
        idx <- as.logical(strand(newgr)=='-')
        GenomicRanges::end(  newgr)[idx] <- GenomicRanges::end(gr)[idx]   - start  # do not switch these lines
        GenomicRanges::start(newgr)[idx] <- GenomicRanges::end(gr)[idx]   - end    # to avoid integrity errors
        formatstr <- '\t\t%d up  flanks: [seqstart%s%d, seqstart%s%d]'
    }
    txt <- sprintf(formatstr, 
                   length(newgr), csign(start), abs(start), csign(end), abs(end))
    
    # Add seq
    if ('seq' %in% names(mcols(gr))){
        assert_is_all_of(bsgenome, 'BSgenome')
        newgr %<>% add_seq(bsgenome)
    }

    # Plot, Message, Return
    if (plot){
        gr$set    <- 'original'
        newgr$set <- 'up flanks'
        allgr <- c(gr, newgr)
        allgr$set %<>% factor(c('original', 'up flanks'))
        plot_intervals(allgr, color_var = 'set', ..., title = txt)
    }
    if (verbose) message(txt)
    newgr
}


#' @rdname up_flank
#' @export
down_flank <- function(
    gr,
    start = 1, 
    end   = 200,
    stranded   = TRUE,
    bsgenome   = NULL,
    verbose    = FALSE,
    plot       = FALSE,
    ...
){
    # Assert
    assert_is_any_of(gr, 'GRanges')
    assert_is_a_number(start)
    assert_is_a_number(end)
    assert_is_a_bool(verbose)
    
    # Flank
    newgr <- gr
    GenomicRanges::start(newgr) <- GenomicRanges::end(newgr) + start
    GenomicRanges::end(newgr)   <- GenomicRanges::end(newgr) + end
    formatstr <- '\t\t%d right flanks : [end%s%d, end%s%d]'
    if (stranded){
      idx <- as.logical(strand(newgr)=='-')
      GenomicRanges::start(newgr)[idx] <- GenomicRanges::start(gr)[idx] - end
      GenomicRanges::end(  newgr)[idx] <- GenomicRanges::start(gr)[idx] - start
      formatstr <- '\t\t%d Down  flanks: [seqstart%s%d, seqstart%s%d]'
    }
    txt <- sprintf(formatstr, length(newgr), csign(start), abs(start), 
                    csign(end), abs(end))
    
    # Add seq
    if ('seq' %in% names(mcols(gr))){
        assert_is_all_of(bsgenome, 'BSgenome')
        newgr %<>% add_seq(bsgenome)
    }
    
    # Plot, Message, Return
    if (plot){
        gr$set <- 'original'
        newgr$set <- 'down flanks'
        allgr <- c(gr, newgr)
        allgr$set %<>% factor(c('original', 'down flanks'))
        plot_intervals(allgr, color_var = 'set', ..., title=txt)
    }
    if (verbose) message(txt)
    newgr
}


#' @rdname up_flank
#' @export
extend <- function(
    gr, 
    start = -22, 
    end  =  22,
    stranded  = TRUE,
    bsgenome  = NULL,
    verbose   = FALSE,
    plot      = FALSE,
    ...
){

    # Assert
    assert_is_any_of(gr, 'GRanges')
    assert_is_a_number(start)
    assert_is_a_number(end)
    assert_is_a_bool(verbose)
    
    # extend
    newgr <- gr
    GenomicRanges::start(newgr) <- GenomicRanges::start(newgr) + start
    GenomicRanges::end(  newgr) <- GenomicRanges::end(  newgr) + end
    formatstr <- '\t\t%d extended ranges: [start%s%d, end%s%d]'
    if (stranded){
      idx <- as.logical(strand(newgr)=='-')
      GenomicRanges::start(newgr)[idx] <- GenomicRanges::start(gr)[idx] - end
      GenomicRanges::end(  newgr)[idx] <- GenomicRanges::end(gr)[idx]   - start
    }
    txt <- sprintf(formatstr, length(newgr), csign(start), abs(start), 
                    csign(end), abs(end))
    
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
        plot_intervals(allgr, color_var = 'set', ..., title=txt)
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
#' # Load
#'     require(magrittr)
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#'     bsinfo <- BSgenome::seqinfo(bsgenome)
#'     
#' # PRNP snp: Kuru resistance variant (G -> T)
#'     gr <- GenomicRanges::GRanges(
#'             'chr20:4699500', strand = '+', seqinfo = bsinfo)
#'     gr %<>% add_seq(bsgenome)
#'     gr %>%  add_inverse_strand(plot = TRUE)
#'     
#' # HBB snp: sickle cell variant (T -> A)
#'     gr <- GenomicRanges::GRanges(
#'             'chr11:5227002-5227002', strand = '-', seqinfo = bsinfo)
#'     gr %<>% add_seq(bsgenome)
#'     gr %>%  add_inverse_strand(plot = TRUE)
#'     
#' # HEXA TATC duplication: Tay-Sachs variant
#'     gr <- GenomicRanges::GRanges(
#'             'chr15:72346580-72346583', strand = '-', seqinfo = bsinfo)
#'     gr %<>% add_seq(bsgenome)
#'     gr %>%  add_inverse_strand(plot = TRUE)
#' @export
add_inverse_strand <- function(gr, verbose = FALSE, plot = FALSE, ...){
    
    # Invert
    complements <- invertStrand(gr)
    
    # Add seq
    if ('seq' %in% names(mcols(gr))){
        complements$seq <- as.character(complement(DNAStringSet(gr$seq)))
    }
    
    # Concatenate
    newgr <- c(gr, complements)
    txt <- sprintf('\t\t%d ranges after adding inverse strands', length(newgr))
    
    # Sort
    newgr <- sortSeqlevels(newgr)
    newgr <- GenomicRanges::sort(newgr)
    
    # Plot
    if (plot){
        gr$set    <- 'sites'
        newgr$set <- 'inv'
        plot_intervals(c(gr, newgr), color_var = 'set', ..., title = txt)
    }
    
    # Message
    if (verbose) cmessage(txt)
    newgr
}


