
summarize_loci <- function(gr){
    sprintf('%s:%s-%s', 
            as.character(GenomeInfoDb::seqnames(gr)), 
            start(gr), 
            end(gr))
}


#' @rdname straddle
#' @export
left_flank <- function(
    gr, 
    leftstart  = -200,
    leftend    = -1,
    plot       = TRUE,
    verbose    = TRUE
){
    # Assert
    assert_is_any_of(gr, 'GRanges')
    assert_is_a_number(leftstart)
    assert_is_a_number(leftend)
    assert_is_a_bool(verbose)
    
    # Flank
    newranges <- gr
    end(newranges)   <- start(gr) + leftend
    start(newranges) <- start(gr) + leftstart
    txt <- sprintf('\t\t%d left  flanks: [start%s%d, start%s%d]', 
                   length(newranges),
                   csign(leftstart), 
                   abs(leftstart), 
                   csign(leftend),
                   abs(leftend))

    # Plot, Message, Return
    if (plot){
        grlist <- GRangesList(sites = gr, leftflanks = newranges)
        plot_intervals(grlist, title = txt)
    }
    if (verbose) message(txt)
    newranges
}


#' @rdname straddle
#' @export
right_flank <- function(
    gr,
    rightstart = 1, 
    rightend   = 200,
    plot       = TRUE,
    verbose    = TRUE
){
    # Assert
    assert_is_any_of(gr, 'GRanges')
    assert_is_a_number(rightstart)
    assert_is_a_number(rightend)
    assert_is_a_bool(verbose)
    
    # Flank
    newranges <- gr
    start(newranges) <- end(newranges) + rightstart
    end(newranges)   <- end(newranges) + rightend
    txt <- sprintf('\t\t%d right flanks : [end%s%d, end%s%d]', 
                    length(newranges),
                    csign(rightstart), 
                    abs(rightstart), 
                    csign(rightend), 
                    abs(rightend))
    
    # Plot, Message, Return
    if (plot){
        grlist <- GRangesList(sites = gr, rightflanks = newranges)
        plot_intervals(grlist, title = txt)
    }
    if (verbose) message(txt)
    newranges
}


#' @rdname straddle
#' @export
extend <- function(
    gr, 
    leftstart = -22, 
    rightend  =  22,
    plot      = TRUE,
    verbose   = TRUE
){

    # Assert
    assert_is_any_of(gr, 'GRanges')
    assert_is_a_number(leftstart)
    assert_is_a_number(rightend)
    assert_is_a_bool(verbose)
    
    # Extend
    newranges <- gr
    start(newranges) <- start(newranges) + leftstart
    end(newranges)   <- end(newranges)   + rightend
    txt <- sprintf('\t\t%d extended ranges: [start%s%d, end%s%d]', 
                    length(newranges),
                    csign(leftstart), 
                    abs(leftstart), 
                    csign(rightend), 
                    abs(rightend))
    
    # Plot, Message, Return
    if (plot){
        grlist <- GRangesList(original = gr, extended = newranges)
        plot_intervals(grlist, title = txt)
    }
    if (verbose) message(txt)
    newranges
    
}


#' Extend or flank GRanges
#' @param gr   \code{\link[GenomicRanges]{GRanges-class}}
#' @param leftstart  number: left start (relative to range start)
#' @param leftend    number: left end   (relative to range start)
#' @param rightstart number: right start (relative to range end)
#' @param rightend   number: right end   (relative to range end)
#' @param plot       logical(1)
#' @param verbose    logical(1)
#' @return a \code{\link[GenomicRanges]{GRanges-class}}
#' @examples 
#' # Read ranges
#'     bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'     gr <- read_bed(bedfile, 'mm10', plot = FALSE)
#' 
#' # Flank/Extend
#'     left_flank( gr, -200,  -1)
#'     right_flank(gr,    1, 200)
#'     extend(     gr, -200, 200)
#'     
#' # Single-verb alternative
#'     straddle(gr, leftstart=-200, leftend=-1)    # left flank
#'     straddle(gr, rightstart=1,   rightend=200)  # right flank
#'     straddle(gr, leftstart=-200, rightend=200)  # extend
#' @export
straddle <- function(
    gr, 
    leftstart  = NULL, 
    leftend    = NULL, 
    rightstart = NULL, 
    rightend   = NULL,
    plot       = TRUE,
    verbose    = TRUE
){
    
    # Extend    
    newranges <- gr
    if (is.numeric(leftstart) & is.numeric(rightend)){
        extend( gr, 
                leftstart  = leftstart, 
                rightend   = rightend,
                plot       = plot, 
                verbose    = verbose)
        if (!is.null(leftend) | !is.null(rightstart)){
            warning('Ignore leftend/rightstart to resolve ambiguity')
        }
    }
    
    # Left flank
    if (is.numeric(leftstart) & is.numeric(leftend)){
        left_flank( gr, 
                    leftstart  = leftstart, 
                    leftend    = leftend, 
                    plot       = plot, 
                    verbose    = verbose)
        if (!is.null(rightstart) | !is.null(rightend)){
            warning('Ignore rightstart/rightend to resolve ambiguity')
        }
    }
    
    # Right flank
    if (is.numeric(rightstart) & is.numeric(rightend)){
        right_flank(gr, 
                    rightstart = rightstart, 
                    rightend   = rightend, 
                    plot       = plot, 
                    verbose    = verbose)
        if (!is.null(leftstart) | !is.null(leftend)){
            warning('Ignore leftstart/leftend to resolve ambiguity')
        }
    }
}


#' Double flank
#' 
#' Flank left and right and merge overlaps
#' @param gr          \code{\link[GenomicRanges]{GRanges-class}}
#' @param leftstart   number: left flank start  (relative to range start)
#' @param leftend     number: left flank  end   (relative to range start)
#' @param rightstart  number: right flank start (relative to range end)
#' @param rightend    number: right flank end   (relative to range end)
#' @param complement  logical(1): whether to add complementary strand
#' @param plot        logical(1)
#' @param verbose     logical(1)
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @examples
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' gr <- read_bed(bedfile, 'mm10', plot = FALSE)
#' double_flank(gr)
#' @export
double_flank <- function(
    gr,
    leftstart  = -200,
    leftend    =   -1,
    rightstart =    1,
    rightend   =  200,
    complement = TRUE,
    plot       = TRUE,
    verbose    = TRUE
){
    # Comply
    . <- NULL
    
    # Flank
    if (verbose) cmessage('\tFlank fourways')
    left <-  left_flank(gr, leftstart,   leftend,  plot=FALSE, verbose=verbose)
    right <- right_flank(gr, rightstart, rightend, plot=FALSE, verbose=verbose)
    newranges <- c(left, right)
    if (verbose) cmessage('\t\t%d combined (left + right)', length(newranges))

    # Plot
    if (plot)  plot_intervals(GRangesList(sites = gr, flanks = newranges))

    # Merge overlaps
    newranges %<>% reduce() # GenomicRanges::reduce
    if (verbose) cmessage('\t\t%d after merging overlaps', length(newranges))
    
    # Return
    newranges
}

