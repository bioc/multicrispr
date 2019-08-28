

#' Left flank 
#' @param granges   GenomicRanges::GRanges
#' @param leftstart left flank start (relative to range start)
#' @param leftend   left flank end   (relative to range start)
#' @param verbose   logical(1)
#' @return GenomicRanges::GRanges
#' @export
#' @examples 
#' require(magrittr)
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' granges <- read_bed(bedfile, bsgenome)
#' granges %>% head(3)
#' granges %>% head(3) %>% left_flank()
left_flank <- function(
    granges, 
    leftstart = -200,
    leftend   = -1,
    verbose   = TRUE
){
    # Assert
    assertive.base::assert_is_identical_to_true(
        methods::is(granges, 'GRanges'))
    assertive.types::assert_is_a_number(leftstart)
    assertive.types::assert_is_a_number(leftend)
    assertive.types::assert_is_a_bool(verbose)
    
    # Flank
    GenomicRanges::end(granges)   <- GenomicRanges::start(granges) + leftend
    GenomicRanges::start(granges) <- GenomicRanges::start(granges) + leftstart
    
    # TODO: Check that they remain within-range!

    # Return
    cmessage('\t\t%d left  flanks : [start%s%d, start%s%d]', 
            length(granges),
            csign(leftstart), 
            abs(leftstart), 
            csign(leftend),
            abs(leftend))

    granges
}


#' Right flank 
#' @param granges      data.table(chr, start, end)
#' @param rightstart flank start relative to range start
#' @param rightend   flank end   relative to range start
#' @param verbose     logical(1)
#' @return GenomicRanges::GRanges
#' @export
#' @examples 
#' require(magrittr)
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' granges <- read_bed(bedfile, bsgenome)
#' granges %>% head(1)
#' granges %>% head(1) %>% right_flank()
#' @export
right_flank <- function(
    granges,
    rightstart = 1, 
    rightend   = 200, 
    verbose     = TRUE
){
    # Assert
    assertive.base::assert_is_identical_to_true(
        methods::is(granges, 'GRanges'))
    assertive.types::assert_is_a_number(rightstart)
    assertive.types::assert_is_a_number(rightend)
    assertive.types::assert_is_a_bool(verbose)
    
    # Flank
    GenomicRanges::start(granges) <- GenomicRanges::end(granges) + rightstart
    GenomicRanges::end(granges)   <- GenomicRanges::end(granges) + rightend
    
    # TODO: check consistency

    # Return
    cmessage('\t\t%d right flanks : [end%s%d, end%s%d]', 
            length(granges),
            csign(rightstart), 
            abs(rightstart), 
            csign(rightend), 
            abs(rightend))

    granges
}


#' Slop (i.e. extend left/right)
#' @param granges   GenomicRanges::GRanges
#' @param leftstart flank start relative to range start
#' @param rightend  flank end   relative to range start
#' @param verbose   logical(1)
#' @return GenomicRanges::GRanges
#' @export
#' @examples 
#' require(magrittr)
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' granges <- read_bed(bedfile, bsgenome)
#' granges %>% head(1)
#' granges %>% head(1) %>% slop()
#' @export
slop <- function(
    granges, 
    leftstart = -22, 
    rightend   =  22, 
    verbose     = TRUE
){

    # Assert
    assertive.base::assert_is_identical_to_true(methods::is(granges, 'GRanges'))
    assertive.types::assert_is_a_number(leftstart)
    assertive.types::assert_is_a_number(rightend)
    assertive.types::assert_is_a_bool(verbose)
    
    # Slop
    GenomicRanges::start(granges) <- GenomicRanges::start(granges) + leftstart
    GenomicRanges::end(granges)   <- GenomicRanges::end(granges)   + rightend
    
    # TODO check consistency

    # Return
    if (verbose) cmessage(
                    '\t\t%d slopped granges: [start%s%d, end%s%d]', 
                    length(granges),
                    csign(leftstart), 
                    abs(leftstart), 
                    csign(rightend), 
                    abs(rightend))
    granges
    
}


#' Complement granges
#' @param granges GenomicRanges::GRanges
#' @param verbose logical(1)
#' @return GenomicRanges::GRanges
#' @examples
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr') 
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' granges <- read_bed(bedfile, bsgenome)
#' granges
#' complement(granges) 
#' @export
complement <- function(granges, verbose = TRUE){

    # Assert
    assertive.base::assert_is_identical_to_true(methods::is(granges, 'GRanges'))
    assertive.sets::assert_is_subset(
        unique(S4Vectors::runValue(GenomicRanges::strand(granges))), 
        c('-', '+'))

    # Complement
    complranges <- granges
    GenomicRanges::strand(complranges) %<>% (function(y) ifelse(y== '+', '-', '+'))

    # Return    
    if (verbose) cmessage('\t\t%d strand-complementary granges', 
                            length(complranges))
    complranges

}


#' Flank fourways
#' 
#' Flank left and right, for both strands, and merge overlaps
#' @param granges     GenomicRanges::GRanges
#' @param leftstart   numeric(1): left flank start  (from range start)
#' @param leftend     numeric(1): left flank  end   (from range start)
#' @param rightstart  numeric(1): right flank start (from range end)
#' @param rightend    numeric(1): right flank end   (from range end)
#' @param verbose     logical(1): report?
#' @return data.table(chr, start, end, strand)
#' @examples
#' require(magrittr)
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' granges <- read_bed(bedfile, bsgenome)
#' granges
#' granges %>% flank_fourways()
##' @export
flank_fourways <- function(
    granges,
    leftstart  = -200,
    leftend    =   -1,
    rightstart =    1,
    rightend   =  200,
    verbose    = TRUE
){
    # Comply
    . <- NULL
    
    # Flank
    if (verbose) cmessage('\tFlank fourways')
    newranges <- c( left_flank(  granges,
                                leftstart  = leftstart, 
                                leftend    = leftend, 
                                verbose    = verbose),
                    right_flank(granges,
                                rightstart = rightstart,
                                rightend   = rightend,
                                verbose    = verbose))
    if (verbose) cmessage('\t\t%d combined (left + right)', 
                            length(newranges))

    # Complement
    newranges %<>% c(complement(., verbose = FALSE))
    if (verbose) cmessage('\t\t%d after adding strand-complements', 
                            length(newranges))

    # Merge overlaps
    newranges %<>% GenomicRanges::reduce()
    if (verbose) cmessage('\t\t%d after merging overlaps', 
                            length(newranges))
    
    # Return
    newranges
}


#' Slop granges for both strands, merging overlaps
#' @param granges   GenomicRanges::GRanges
#' @param leftstart numeric(1)
#' @param rightend  numeric(1)
#' @param verbose   logical(1)
#' @return GenomicRanges::GRanges
#' @examples
#' require(magrittr)
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' granges <- read_bed(bedfile, bsgenome)
#' granges
#' granges %>% slop_fourways()
#' @export
slop_fourways <- function(
    granges,
    leftstart = -22,
    rightend   = 22,
    verbose     = TRUE
){
    # Comply
    . <- NULL
    
    # Slop
    if (verbose) cmessage('\tSlop fourways')
    newranges  <-   granges %>% 
                    slop(   leftstart = leftstart, 
                            rightend  = rightend, 
                            verbose   = verbose)
    
    # Complement
    newranges %<>% c(complement(., verbose = FALSE))
    if (verbose)   cmessage('\t\t%d after adding strand-complements', 
                            length(newranges))
    
    # Merge overlaps
    newranges %<>% GenomicRanges::reduce()
    if (verbose) cmessage('\t\t%d after merging overlaps', 
                            length(newranges))
    
    # Return
    newranges
}


