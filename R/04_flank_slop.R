
summarize_loci <- function(gr){
    sprintf('%s:%s-%s', 
            as.character(seqnames(gr)), 
            start(gr), 
            end(gr))
}


#' Left flank 
#' @param gr   \code{\link[GenomicRanges]{GRanges-class}}
#' @param leftstart number: flank start (relative to range start)
#' @param leftend   number: flank end   (relative to range start)
#' @param plot      logical(1)
#' @param verbose   logical(1)
#' @return a \code{\link[GenomicRanges]{GRanges-class}}
#' @export
#' @examples 
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' gr <- read_bed(bedfile, bsgenome, plot = FALSE)
#' left_flank(gr)
left_flank <- function(
    gr, 
    leftstart = -200,
    leftend   = -1,
    plot      = TRUE,
    verbose   = TRUE
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

    # Plot
    if (plot)  plot_intervals(GRangesList(sites = gr, leftflanks = newranges))
    
    # Return
    cmessage('\t\t%d left  flanks : [start%s%d, start%s%d]', 
            length(newranges),
            csign(leftstart), 
            abs(leftstart), 
            csign(leftend),
            abs(leftend))

    newranges
}


#' Right flank 
#' @param gr    \code{\link[GenomicRanges]{GRanges-class}}
#' @param rightstart number: flank start (relative to range end)
#' @param rightend   number: flank end   (relative to range end)
#' @param plot       logical(1)
#' @param verbose    logical(1)
#' @return     \code{\link[GenomicRanges]{GRanges-class}}
#' @export
#' @examples 
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' gr <- read_bed(bedfile, bsgenome, plot = FALSE)
#' right_flank(gr)
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
    
    # Plot
    if (plot)  plot_intervals(GRangesList(sites = gr, rightflanks = newranges))
    
    # Return
    cmessage('\t\t%d right flanks : [end%s%d, end%s%d]', 
            length(newranges),
            csign(rightstart), 
            abs(rightstart), 
            csign(rightend), 
            abs(rightend))
    newranges
}


#' Slop (i.e. extend left/right)
#' @param gr        \code{\link[GenomicRanges]{GRanges-class}}
#' @param leftstart number: flank start (relative to range start)
#' @param rightend  number: flank end   (relative to range end)
#' @param plot      logical(1)
#' @param verbose   logical(1)
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @export
#' @examples 
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' gr <- read_bed(bedfile, bsgenome, plot = FALSE)
#' slop(gr)
#' @export
slop <- function(
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
    
    # Slop
    newranges <- gr
    start(newranges) <- start(newranges) + leftstart
    end(newranges)   <- end(newranges)   + rightend
    
    # Plot
    if (plot)  plot_intervals(GRangesList(original = gr, slopped = newranges))

    
    # Return
    if (verbose) cmessage(
                    '\t\t%d slopped ranges: [start%s%d, end%s%d]', 
                    length(newranges),
                    csign(leftstart), 
                    abs(leftstart), 
                    csign(rightend), 
                    abs(rightend))
    newranges
    
}


#' Flank fourways
#' 
#' Flank left and right, for both strands, and merge overlaps
#' @param gr          \code{\link[GenomicRanges]{GRanges-class}}
#' @param leftstart   number: left flank start  (relative to range start)
#' @param leftend     number: left flank  end   (relative to range start)
#' @param rightstart  number: right flank start (relative to range end)
#' @param rightend    number: right flank end   (relative to range end)
#' @param plot        logical(1)
#' @param verbose     logical(1)
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @examples
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' gr <- read_bed(bedfile, bsgenome, plot = FALSE)
#' flank_fourways(gr)
#' @export
flank_fourways <- function(
    gr,
    leftstart  = -200,
    leftend    =   -1,
    rightstart =    1,
    rightend   =  200,
    plot       = TRUE,
    verbose    = TRUE
){
    # Comply
    . <- NULL
    
    # Flank
    if (verbose) cmessage('\tFlank fourways')
    left <-  left_flank( gr, leftstart, leftend, plot = FALSE, verbose=verbose)
    right <- right_flank(gr,rightstart, rightend, plot= FALSE, verbose=verbose)
    newranges <- c(left, right)
    if (verbose) cmessage('\t\t%d combined (left + right)', length(newranges))

    # Complement
    newranges %<>% c(invertStrand(.))
    if (verbose) cmessage('\t\t%d after adding strand inverts',
                            length(newranges))

    # Plot
    if (plot)  plot_intervals(GRangesList(sites = gr, flanks = newranges))

    # Merge overlaps
    newranges %<>% reduce() # GenomicRanges::reduce
    if (verbose) cmessage('\t\t%d after merging overlaps', length(newranges))
    
    # Return
    newranges
}


#' Slop granges for both strands, merging overlaps
#' @param gr   \code{\link[GenomicRanges]{GRanges-class}}
#' @param leftstart number
#' @param rightend  number
#' @param plot      logical(1)
#' @param verbose   logical(1)
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @examples
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' gr <- read_bed(bedfile, bsgenome, plot = FALSE)
#' slop_fourways(gr)
#' @export
slop_fourways <- function(
    gr,
    leftstart = -22,
    rightend  =  22,
    plot      = TRUE, 
    verbose   = TRUE
){
    # Comply
    . <- NULL
    
    # Slop
    if (verbose) cmessage('\tSlop fourways')
    newranges <- slop(gr, leftstart, rightend, verbose = verbose)
    
    # Complement
    newranges %<>% c(invertStrand(.))
    if (verbose)   cmessage('\t\t%d after adding strand inverts', 
                            length(newranges))
    
    # Merge overlaps
    newranges %<>% reduce()
    if (verbose) cmessage('\t\t%d after merging overlaps', length(newranges))

    # Plot
    if (plot)  plot_intervals(GRangesList(sites = gr, slopped = newranges))

    # Return
    newranges
}

