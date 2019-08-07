

#' Left flank 
#' @param ranges      data.table(chr, start, end, strand)
#' @param bsgenome    BSgenome, e.g. BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' @param leftstart   left flank start (relative to range start)
#' @param leftend     left flank end   (relative to range start)
#' @param verbose     logical(1)
#' @return data.table
#' @export
#' @examples 
#' require(magrittr)
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#' ranges <- read_bed(bedfile)
#' ranges %>% head(3)
#' ranges %>% head(3) %>% left_flank(bsgenome)
left_flank <- function(
    ranges, 
    bsgenome,
    leftstart = -200,
    leftend   = -1,
    verbose   = TRUE
){
    # Assert
    assertive.types::assert_is_data.table(ranges)
    assertive.sets::assert_is_subset(c('chr', 'start', 'end'), names(ranges))
    assertive.types::assert_is_a_number(leftstart)
    assertive.types::assert_is_a_number(leftend)
    assertive.types::assert_is_a_bool(verbose)
    assertive.base::assert_is_identical_to_true(
        methods::is(bsgenome, 'BSgenome'))
    
    # Flank
    chr <- start <- end <- chrlength <- NULL
    newranges <-  data.table::copy(ranges) %>% 
                    extract(, end   := start + leftend)   %>%  # dont switch
                    extract(, start := start + leftstart) %>%  # lines!
                    extract()
    newranges [ , chrlength  := GenomeInfoDb::seqlengths(bsgenome)[chr] ]
    tmp <- newranges [, assertive.base::assert_all_are_true(start >= 1)      ]
    tmp <- newranges [, assertive.base::assert_all_are_true(end <= chrlength)]
    newranges[, 'chrlength' := NULL]

    # Return
    cmessage('\t\t%d left  flanks : [start%s%d, start%s%d]', 
            nrow(newranges),
            csign(leftstart), 
            abs(leftstart), 
            csign(leftend),
            abs(leftend))

    newranges[]
}


#' Right flank 
#' @param ranges      data.table(chr, start, end)
#' @param bsgenome    BSgenome, e.g. BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' @param rightstart flank start relative to range start
#' @param leftend   flank end   relative to range start
#' @param verbose     logical(1)
#' @return data.table
#' @export
#' @examples 
#' require(magrittr)
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#' ranges <- read_bed(bedfile)
#' ranges %>% head(1)
#' ranges %>% head(1) %>% right_flank(bsgenome)
#' @export
right_flank <- function(
    ranges,
    bsgenome,
    rightstart = 1, 
    rightend   = 200, 
    verbose     = TRUE
){
    # Assert
    assertive.types::assert_is_data.table(ranges)
    assertive.sets::assert_is_subset(c('chr', 'start', 'end'), names(ranges))
    assertive.types::assert_is_a_number(rightstart)
    assertive.types::assert_is_a_number(rightend)
    assertive.types::assert_is_a_bool(verbose)
    
    # Flank
    chr <- start <- end <- chrlength <- NULL
    newranges  <-   data.table::copy(ranges)               %>% 
                    extract(, start := end + rightstart)  %>% 
                    extract(, end   := end + rightend)    %>% 
                    extract()
    newranges [ , chrlength  := GenomeInfoDb::seqlengths(bsgenome)[chr] ]
    tmp <- newranges [, assertive.base::assert_all_are_true(start >= 1)      ]
    tmp <- newranges [, assertive.base::assert_all_are_true(end <= chrlength)]
    newranges[, 'chrlength' := NULL]

    # Return
    cmessage('\t\t%d right flanks : [end%s%d, end%s%d]', 
            nrow(newranges),
            csign(rightstart), 
            abs(rightstart), 
            csign(rightend), 
            abs(rightend))

    newranges[]
}


#' Slop (i.e. extend left/right)
#' @param ranges      data.table(chr, start, end)
#' @param bsgenome    BSgenome, e.g. BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' @param leftstart flank start relative to range start
#' @param rightend   flank end   relative to range start
#' @param verbose     logical(1)
#' @return data.table
#' @export
#' @examples 
#' require(magrittr)
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#' ranges <- read_bed(bedfile)
#' ranges %>% head(1)
#' ranges %>% head(1) %>% slop(bsgenome)
#' @export
slop <- function(
    ranges, 
    bsgenome,
    leftstart = -22, 
    rightend   =  22, 
    verbose     = TRUE
){

    # Assert
    assertive.types::assert_is_data.table(ranges)
    assertive.sets::assert_is_subset(c('chr', 'start', 'end'), names(ranges))
    assertive.types::assert_is_a_number(leftstart)
    assertive.types::assert_is_a_number(rightend)
    assertive.types::assert_is_a_bool(verbose)
    
    # Slop
    chr <- start <- end <- chrlength <- NULL
    newranges  <-   data.table::copy(ranges)                 %>% 
                    extract(, start := start + leftstart)  %>%  
                    extract(, end   := end   + rightend)    %>% 
                    extract()
    newranges [ , chrlength  := GenomeInfoDb::seqlengths(bsgenome)[chr] ]
    tmp <- newranges [, assertive.base::assert_all_are_true(start >= 1)      ]
    tmp <- newranges [, assertive.base::assert_all_are_true(end <= chrlength)]
    newranges[, 'chrlength' := NULL]

    # Return
    if (verbose) cmessage(
                    '\t\t%d slopped ranges: [start%s%d, end%s%d]', 
                    nrow(newranges),
                    csign(leftstart), 
                    abs(leftstart), 
                    csign(rightend), 
                    abs(rightend))
    newranges[]
    
}

#' Merge overlaps
#' @param x   data.table(chr, start, end, strand)
#' @param bsgenome BSgenome, e.g. BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' @param verbose  logical(1)
#' @return data.table(chr, start, end)
#' @examples 
#' require(magrittr)
#' (rangesdt <- data.table::data.table(
#'             chr    = rep('chr1', 4), 
#'             start  = c(1,    5,   5,  20), 
#'             end    = c(10,  15,  15,  30), 
#'             strand = c('+', '+', '-', '-')))
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' rangesdt %>% merge_overlaps(bsgenome)
#' @export
merge_overlaps <- function(x, bsgenome, verbose = TRUE){
    reduced_ranges <-   as.granges(x, bsgenome) %>% 
                        GenomicRanges::reduce() %>% 
                        as.data.table()
    if (verbose)  cmessage('\t\t%d ranges after merging overlaps', 
                            nrow(reduced_ranges))
    reduced_ranges[]

# Alternative data.table solution:    
# https://stackoverflow.com/a/41748171
}


#' Complement ranges
#' @param ranges data.table(chr, start, end, strand)
#' @param verbose logical(1)
#' @return data.table(chr, start, end, strand)
#' @examples
#' bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex') 
#' ranges <- read_bed(bedfile)
#' complement(ranges) 
#' @export
complement <- function(ranges, verbose = TRUE){

    # Assert
    assertive.types::assert_is_data.table(ranges)
    assertive.sets::assert_is_subset(
        c('chr', 'start', 'end', 'strand'), names(ranges))
    
    # Complement
    strand <- NULL
    assertive.sets::assert_is_subset(unique(ranges$strand), c('-', '+'))

    newranges <-    data.table::copy(ranges) %>% 
                    extract(, strand := ifelse(strand=='+', '-', '+')) %>% 
                    extract()
    
    # Return
    if (verbose) cmessage('\t\t%d strand-complementary ranges', 
                            nrow(newranges))
    newranges[]

}


#' Flank fourways
#' 
#' Flank left and right, for both strands, and merge overlaps
#' @param ranges      data.table(chr, start, end, strand)
#' @param bsgenome    BSgenome, e.g. BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' @param leftstart   numeric(1): left flank start  (from range start)
#' @param leftend     numeric(1): left flank  end   (from range start)
#' @param rightstart  numeric(1): right flank start (from range end)
#' @param rightend    numeric(1): right flank end   (from range end)
#' @param verbose     logical(1): report?
#' @return data.table(chr, start, end, strand)
#' @examples 
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#' ranges <- read_bed(bedfile)
#' flank_fourways(ranges, bsgenome)
##' @export
flank_fourways <- function(
    ranges,
    bsgenome,
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
    newranges <- rbind(left_flank(
                            ranges,
                            bsgenome   = bsgenome, 
                            leftstart  = leftstart, 
                            leftend    = leftend, 
                            verbose    = verbose
                        ),
                        right_flank(
                            ranges,
                            bsgenome   = bsgenome,
                            rightstart = rightstart,
                            rightend   = rightend,
                            verbose    = verbose
                        )
                )
    if (verbose) cmessage('\t\t%d ranges combined (left + right)', 
                            nrow(newranges))

    # Complement
    newranges %<>% rbind(complement(., verbose = FALSE))
    newranges %>% data.table::setorderv(c('chr', 'start', 'end'))
    if (verbose) cmessage('\t\t%d ranges after adding strand-complements', 
                            nrow(newranges))

    # Merge overlaps
    newranges %>% merge_overlaps(bsgenome = bsgenome, verbose = verbose) 
}


#' Slop ranges for both strands, merging overlaps
#' @param ranges      data.table(chr, start, end, strand)
#' @param bsgenome    BSgenome, e.g. BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' @param leftstart numeric(1)
#' @param rightend   numeric(1)
#' @param verbose     logical(1)
#' @return data.table(chr, start, end, strand)
#' @examples
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#' ranges <- read_bed(bedfile)
#' slop_fourways(ranges, bsgenome)
#' @export
slop_fourways <- function(
    ranges,
    bsgenome,
    leftstart = -22, 
    rightend   = 22, 
    verbose     = TRUE
){
    # Comply
    . <- NULL
    
    # Slop
    if (verbose) cmessage('\tSlop fourways')
    newranges <-    ranges %>% 
                    slop(
                        bsgenome  = bsgenome,
                        leftstart = leftstart, 
                        rightend  = rightend, 
                        verbose   = verbose
                    )
    
    # Complement
    newranges %<>% rbind(complement(., verbose = FALSE))
    if (verbose)   cmessage('\t\t%d ranges after adding strand-complements', 
                            nrow(newranges))
    
    # Merge overlaps
    newranges %>% merge_overlaps(bsgenome = bsgenome, verbose = verbose)
}


