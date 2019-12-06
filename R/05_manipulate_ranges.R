#' Add sequence to GRanges
#' @param gr        \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome  \code{\link[BSgenome]{BSgenome-class}}
#' @param verbose   logical(1)
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
add_seq <- function(gr, bsgenome, verbose = TRUE){
    
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
                        as.character = TRUE))
    
    # Return
    gr
}


summarize_loci <- function(gr){
    sprintf('%s:%s-%s', 
            as.character(seqnames(gr)), 
            start(gr), 
            end(gr))
}


#' Flank or Extend GRanges
#' @param gr   \code{\link[GenomicRanges]{GRanges-class}}
#' @param leftstart  number: left start (relative to range start)
#' @param leftend    number: left end   (relative to range start)
#' @param rightstart number: right start (relative to range end)
#' @param rightend   number: right end   (relative to range end)
#' @param bsgenome   NULL (default) or \code{\link[BSgenome]{BSgenome-class}}.
#'                   Required to update gr$seq if present.
#' @param verbose    TRUE (default) or FALSE
#' @param plot       TRUE (default) or FALSE
#' @param ...        passed to \code{\link{plot_intervals}}
#' @return a \code{\link[GenomicRanges]{GRanges-class}}
#' @examples 
#' # Small GRanges
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#'     gr <- GenomicRanges::GRanges(
#'             c(HBB = 'chr11:5227002-5227002', PRNP = 'chr20:4699500'), 
#'             strand = c('-', '+'), 
#'             seqinfo = BSgenome::seqinfo(bsgenome))
#'     gr %>% left_flank(-22, -1)
#'     gr %>% right_flank( 1, 22)
#'     gr %>%   up_flank(-22, -1)
#'     gr %>% down_flank(1, 22)
#'     gr %>% extend(-22, 22)
#' 
#' # Large GRanges
#'     require(magrittr)
#'     bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'     gr <- bed_to_granges(bedfile, 'mm10', plot = FALSE)
#'     gr %>%   left_flank(-200,  -1)
#'     gr %>%  right_flank(   1, 200)
#'     gr %>%       extend(-200, 200)
#' @seealso \code{\link{straddle}} (single verb function encompassing all of 
#'          left_flank, right_flank, and extend) and \code{\link{double_flank}}.
#' @export
left_flank <- function(
    gr, 
    leftstart  = -200,
    leftend    = -1,
    bsgenome   = NULL,
    verbose    = TRUE,
    plot       = TRUE,
    ...
){
    # Assert
    assert_is_any_of(gr, 'GRanges')
    assert_is_a_number(leftstart)
    assert_is_a_number(leftend)
    assert_is_a_bool(verbose)
    
    # Flank
    newgr <- gr
    start(newgr) <- start(gr) + leftstart
    end(newgr)   <- start(gr) + leftend
    txt <- sprintf('\t\t%d left  flanks: [start%s%d, start%s%d]', 
                    length(newgr), csign(leftstart), abs(leftstart), 
                    csign(leftend), abs(leftend))
    
    # Add seq
    if ('seq' %in% names(mcols(gr))){
        assert_is_all_of(bsgenome, 'BSgenome')
        newgr %<>% add_seq(bsgenome)
    }

    # Plot, Message, Return
    if (plot){
        gr$set    <- 'original'
        newgr$set <- 'leftflanks'
        allgr <- c(gr, newgr)
        allgr$set %<>% factor(c('original', 'leftflanks'))
        plot_intervals(allgr, color_var = 'set', ..., title = txt)
    }
    if (verbose) message(txt)
    newgr
}

#' @rdname left_flank
#' @export
up_flank <- function(
    gr, 
    upstart    = -200,
    upend      = -1,
    bsgenome   = NULL,
    verbose    = TRUE,
    plot       = TRUE,
    ...
){
    # Assert
    assert_is_any_of(gr, 'GRanges')
    assert_is_a_number(upstart)
    assert_is_a_number(upend)
    assert_is_a_bool(verbose)
    
    # Left flank "+" ranges
    newgr <- gr
    idx <- as.logical(strand(newgr)=='+')
    start(newgr)[idx] <- start(gr)[idx] + upstart
      end(newgr)[idx] <- start(gr)[idx] + upend
      
    # Right flank "-"  ranges
    idx <- as.logical(strand(newgr)=='-')
      end(newgr)[idx] <- end(gr)[idx]   - upstart  # do not switch these lines
    start(newgr)[idx] <- end(gr)[idx]   - upend    # to avoid integrity errors
    txt <- sprintf('\t\t%d up  flanks: [seqstart%s%d, seqstart%s%d]', 
                    length(newgr), csign(upstart), abs(upstart), 
                    csign(upend), abs(upend))
    
    # Add seq
    if ('seq' %in% names(mcols(gr))){
        assert_is_all_of(bsgenome, 'BSgenome')
        newgr %<>% add_seq(bsgenome)
    }

    # Plot, Message, Return
    if (plot){
        gr$set    <- 'original'
        newgr$set <- 'upflanks'
        allgr <- c(gr, newgr)
        allgr$set %<>% factor(c('original', 'upflanks'))
        plot_intervals(allgr, color_var = 'set', ..., title = txt)
    }
    if (verbose) message(txt)
    newgr
}



#' @rdname left_flank
#' @export
right_flank <- function(
    gr,
    rightstart = 1, 
    rightend   = 200,
    bsgenome   = NULL,
    verbose    = TRUE,
    plot       = TRUE,
    ...
){
    # Assert
    assert_is_any_of(gr, 'GRanges')
    assert_is_a_number(rightstart)
    assert_is_a_number(rightend)
    assert_is_a_bool(verbose)
    
    # Flank
    newgr <- gr
    start(newgr) <- end(newgr) + rightstart
    end(newgr)   <- end(newgr) + rightend
    txt <- sprintf('\t\t%d right flanks : [end%s%d, end%s%d]', 
                    length(newgr),
                    csign(rightstart), 
                    abs(rightstart), 
                    csign(rightend), 
                    abs(rightend))
    
    # Add seq
    if ('seq' %in% names(mcols(gr))){
        assert_is_all_of(bsgenome, 'BSgenome')
        newgr %<>% add_seq(bsgenome)
    }
    
    # Plot, Message, Return
    if (plot){
        gr$set <- 'original'
        newgr$set <- 'rightflanks'
        allgr <- c(gr, newgr)
        allgr$set %<>% factor(c('original', 'rightflanks'))
        plot_intervals(allgr, color_var = 'set', ..., title=txt)
    }
    if (verbose) message(txt)
    newgr
}

#' @rdname left_flank
#' @export
down_flank <- function(
    gr, 
    downstart    = +1,
    downend      = +200,
    bsgenome   = NULL,
    verbose    = TRUE,
    plot       = TRUE,
    ...
){
    # Assert
    assert_is_any_of(gr, 'GRanges')
    assert_is_a_number(downstart)
    assert_is_a_number(downend)
    assert_is_a_bool(verbose)
    
    # Flank
    newgr <- gr
    idx <- as.logical(strand(newgr)=='+')
    start(newgr)[idx] <- end(gr)[idx] + downstart
      end(newgr)[idx] <- end(gr)[idx] + downend
      
    idx <- as.logical(strand(newgr)=='-')
    start(newgr)[idx] <- start(gr)[idx] - downend
      end(newgr)[idx] <- start(gr)[idx] - downstart
    txt <- sprintf('\t\t%d Down  flanks: [seqstart%s%d, seqstart%s%d]', 
                    length(newgr), csign(downstart), abs(downstart), 
                    csign(downend), abs(downend))
    
    # Add seq
    if ('seq' %in% names(mcols(gr))){
        assert_is_all_of(bsgenome, 'BSgenome')
        newgr %<>% add_seq(bsgenome)
    }

    # Plot, Message, Return
    if (plot){
        gr$set    <- 'original'
        newgr$set <- 'downflanks'
        allgr <- c(gr, newgr)
        allgr$set %<>% factor(c('original', 'downflanks'))
        plot_intervals(allgr, color_var = 'set', ..., title = txt)
    }
    if (verbose) message(txt)
    newgr
}



#' @rdname left_flank
#' @export
extend <- function(
    gr, 
    leftstart = -22, 
    rightend  =  22,
    bsgenome  = NULL,
    verbose   = TRUE,
    plot      = TRUE,
    ...
){

    # Assert
    assert_is_any_of(gr, 'GRanges')
    assert_is_a_number(leftstart)
    assert_is_a_number(rightend)
    assert_is_a_bool(verbose)
    
    # Extend
    newgr <- gr
    start(newgr) <- start(newgr) + leftstart
    end(newgr)   <- end(newgr)   + rightend
    txt <- sprintf('\t\t%d extended ranges: [start%s%d, end%s%d]', 
                    length(newgr),
                    csign(leftstart), 
                    abs(leftstart), 
                    csign(rightend), 
                    abs(rightend))
    
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


#' Straddle GRanges
#' @param gr   \code{\link[GenomicRanges]{GRanges-class}}
#' @param leftstart  number: left start (relative to range start)
#' @param leftend    number: left end   (relative to range start)
#' @param rightstart number: right start (relative to range end)
#' @param rightend   number: right end   (relative to range end)
#' @param bsgenome   \code{\link[BSgenome]{BSgenome-class}}
#'                   Required to update gr$seq if present.
#' @param verbose    TRUE (default) or FALSE
#' @param plot       TRUE (default) or FALSE
#' @param ...         \code{\link{plot_intervals}} arguments
#' @return a \code{\link[GenomicRanges]{GRanges-class}}
#' @examples 
#' # HBB snp: sickle cell variant (T -> A)
#'     require(magrittr)
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#'     bsinfo <- BSgenome::seqinfo(bsgenome)
#'     gr <- GenomicRanges::GRanges(
#'              'chr11:5227002-5227002', strand = '-', seqinfo = bsinfo)
#'     (gr %<>% add_seq(bsgenome))
#'     gr %>% straddle(leftstart = -22, rightend = 22, bsgenome = bsgenome)
#' 
#' # SRF binding sites
#'     bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'     gr <- bed_to_granges(bedfile, 'mm10', plot = FALSE)
#'     straddle(gr, leftstart=-200, leftend=-1)    # left flank
#'     straddle(gr, rightstart=1,   rightend=200)  # right flank
#'     straddle(gr, leftstart=-200, rightend=200)  # extend
#' @author Aditya Bhagwat, after discussing with Michael Lawrence on bioc-devel
#' @export
straddle <- function(
    gr, leftstart  = NULL, leftend = NULL, rightstart = NULL, rightend = NULL, 
    bsgenome = NULL, verbose = TRUE, plot = TRUE, ...
){
    
    # Extend    
    newgr <- gr
    if (is.numeric(leftstart) & is.numeric(rightend)){
        newgr <- extend(
                    gr, leftstart = leftstart, rightend = rightend,
                    bsgenome = bsgenome, verbose = verbose, plot = plot, ...)
        if (!is.null(leftend) | !is.null(rightstart)){
            warning('Ignore leftend/rightstart to resolve ambiguity')
        }
    }
    
    # Left flank
    if (is.numeric(leftstart) & is.numeric(leftend)){
        newgr <- left_flank(
                    gr, leftstart = leftstart, leftend = leftend, 
                    bsgenome = bsgenome, verbose = verbose, plot = plot, ...)
        if (!is.null(rightstart) | !is.null(rightend)){
            warning('Ignore rightstart/rightend to resolve ambiguity')
        }
    }
    
    # Right flank
    if (is.numeric(rightstart) & is.numeric(rightend)){
        newgr <- right_flank(
                    gr, rightstart = rightstart, rightend = rightend, 
                    bsgenome = bsgenome, verbose = verbose, plot = plot, ...)
        if (!is.null(leftstart) | !is.null(leftend)){
            warning('Ignore leftstart/leftend to resolve ambiguity')
        }
    }
    
    # Add seq
    if ('seq' %in% names(mcols(gr))){
        assert_is_all_of(bsgenome, 'BSgenome')
        newgr %<>% add_seq(bsgenome)
    }
    
    # Return
    newgr

}


#' Double flank
#' 
#' Flank left and right and merge overlaps
#' @param gr          \code{\link[GenomicRanges]{GRanges-class}}
#' @param leftstart   number: left flank start  (relative to range start)
#' @param leftend     number: left flank  end   (relative to range start)
#' @param rightstart  number: right flank start (relative to range end)
#' @param rightend    number: right flank end   (relative to range end)
#' @param bsgenome    \code{\link[BSgenome]{BSgenome-class}}
#'                   Required to update gr$seq if present.
#' @param verbose    TRUE (default) or FALSE
#' @param plot       TRUE (default) or FALSE
#' @param ...         \code{\link{plot_intervals}} arguments
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @examples
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' gr <- bed_to_granges(bedfile, 'mm10', plot = FALSE)
#' double_flank(gr)
#' @export
double_flank <- function(
    gr,
    leftstart  = -200,
    leftend    =   -1,
    rightstart =    1,
    rightend   =  200,
    bsgenome   = NULL,
    verbose    = TRUE,
    plot       = TRUE,
    ...
){
    # Comply
    . <- NULL
    
    # Flank
    if (verbose) cmessage('\tFlank fourways')
    left <-  left_flank(gr, leftstart,   leftend,  verbose=verbose, plot=FALSE)
    right <- right_flank(gr, rightstart, rightend, verbose=verbose, plot=FALSE)
    newgr <- c(left, right)
    if (verbose) cmessage('\t\t%d combined (left + right)', length(newgr))

    # Plot
    if (plot){
        gr$set <- 'sites'
        newgr$set <- 'flanks'
        allgr <- c(gr, newgr); allgr$set %<>% factor('original', 'flanks')
        plot_intervals(allgr, color_var = 'set', ...)
    }

    # Merge overlaps
    newgr %<>% GenomicRanges::reduce()
    if (verbose) cmessage('\t\t%d after merging overlaps', length(newgr))
    
    # Add seq
    if ('seq' %in% names(mcols(gr))){
        assert_is_all_of(bsgenome, 'BSgenome')
        newgr %<>% add_seq(bsgenome)
    }
    
    # Return
    newgr
}


#' Add inverse strand
#' @param gr         \code{\link[GenomicRanges]{GRanges-class}}
#' @param verbose    TRUE (default) or FALSE
#' @param plot       TRUE (default) or FALSE
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
#'     gr %>%  add_inverse_strand()
#'     
#' # HBB snp: sickle cell variant (T -> A)
#'     gr <- GenomicRanges::GRanges(
#'             'chr11:5227002-5227002', strand = '-', seqinfo = bsinfo)
#'     gr %<>% add_seq(bsgenome)
#'     gr %>%  add_inverse_strand()
#'     
#' # HEXA TATC duplication: Tay-Sachs variant
#'     gr <- GenomicRanges::GRanges(
#'             'chr15:72346580-72346583', strand = '-', seqinfo = bsinfo)
#'     gr %<>% add_seq(bsgenome)
#'     gr %>%  add_inverse_strand()
#' @export
add_inverse_strand <- function(gr, verbose = TRUE, plot = TRUE, ...){
    
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
