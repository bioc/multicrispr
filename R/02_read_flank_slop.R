#' Read bedfile as data.table
#' @param bedfile        file path
#' @param verbose        logical(1)
#' @param rm_duplicates  logical(1)
#' @return data.table(chr, start, end, strand) 
#' @examples
#' bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#' ranges <- read_bed(bedfile)
#' @importFrom  data.table  :=
#' @export
read_bed <- function(
    bedfile, 
    verbose       = TRUE, 
    rm_duplicates = TRUE
){
    # Assert
    assertive.files::assert_all_are_existing_files(bedfile)
    assertive.types::assert_is_a_bool(verbose)
    assertive.types::assert_is_a_bool(rm_duplicates)

    # Comply
    chr <- start <- end <- strand <- .N <- gap <- width <- NULL

    # Read
    if (verbose) cmessage('\tRead %s', bedfile)
    dt <- data.table::fread(
            bedfile,
            select    = c(seq_len(3), 6),
            col.names = c('chr', 'start', 'end', 'strand'))
    dt %>% data.table::setorderv(c('chr', 'start', 'end', 'strand'))
    if (verbose) cmessage('\t%d ranges on %d chromosomes', 
                            nrow(dt), length(unique(dt$chr)))
    
    # Drop duplicates
    if (rm_duplicates){
        is_duplicated <- cduplicated(dt)
        if (any(is_duplicated)){
            if (verbose) cmessage('\t\t%d after removing duplicates')
            dt %<>% extract(!duplicated)
        }
    }
        
    
    # Report width & gap statistics
    if (verbose){
        dt [ , width := end-start+1]
        dt [ , gap := c(start[2:.N]-end[seq_len(.N-1)], Inf), by = chr ]
        cmessage('\t\t%s NT wide', num2scalarstr(dt$width))
        cmessage('\t\t%s NT apart', num2scalarstr(dt$gap))
        dt [ , c('gap', 'width') := NULL]
    }
    
    # Return
    return(dt[])
}


#' Left flank 
#' @param ranges      data.table(chr, start, end)
#' @param bsgenome    BSgenome, e.g. BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' @param startoffset flank start relative to range start
#' @param endoffset   flank end   relative to range start
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
    startoffset = -200, 
    endoffset   = -1, 
    verbose     = TRUE
){
    # Assert
    assertive.types::assert_is_data.table(ranges)
    assertive.sets::assert_is_subset(c('chr', 'start', 'end'), names(ranges))
    assertive.types::assert_is_a_number(startoffset)
    assertive.types::assert_is_a_number(endoffset)
    assertive.types::assert_is_a_bool(verbose)
    assertive.base::assert_is_identical_to_true(is(bsgenome, 'BSgenome'))
    
    # Flank
    start <- end <- NULL
    newranges <-  data.table::copy(ranges) %>% 
                    extract(, end   := start + endoffset)   %>%  # dont switch
                    extract(, start := start + startoffset) %>%  # lines!
                    extract()
    newranges [ , chrlength  := GenomeInfoDb::seqlengths(bsgenome)[chr] ]
    tmp <- newranges [, assertive.base::assert_all_are_true(start >= 1)      ]
    tmp <- newranges [, assertive.base::assert_all_are_true(end <= chrlength)]
    newranges[, 'chrlength' := NULL]

    # Return
    cmessage('\t\t%d left  flanks : [start%s%d, start%s%d]', 
            nrow(newranges),
            csign(startoffset), 
            abs(startoffset), 
            csign(endoffset),
            abs(endoffset))

    return(newranges)
}


#' Right flank 
#' @param ranges      data.table(chr, start, end)
#' @param bsgenome    BSgenome, e.g. BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' @param startoffset flank start relative to range start
#' @param endoffset   flank end   relative to range start
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
    startoffset = 1, 
    endoffset   = 200, 
    verbose     = TRUE
){
    # Assert
    assertive.types::assert_is_data.table(ranges)
    assertive.sets::assert_is_subset(c('chr', 'start', 'end'), names(ranges))
    assertive.types::assert_is_a_number(startoffset)
    assertive.types::assert_is_a_number(endoffset)
    assertive.types::assert_is_a_bool(verbose)
    
    # Flank
    start <- end <- NULL
    newranges  <-   data.table::copy(ranges)               %>% 
                    extract(, start := end + startoffset)  %>% 
                    extract(, end   := end + endoffset)    %>% 
                    extract()
    newranges [ , chrlength  := GenomeInfoDb::seqlengths(bsgenome)[chr] ]
    tmp <- newranges [, assertive.base::assert_all_are_true(start >= 1)      ]
    tmp <- newranges [, assertive.base::assert_all_are_true(end <= chrlength)]
    newranges[, 'chrlength' := NULL]
   
    # Return
    cmessage('\t\t%d right flanks : [end%s%d, end%s%d]', 
            nrow(newranges),
            csign(startoffset), 
            abs(startoffset), 
            csign(endoffset), 
            abs(endoffset))

    return(newranges)
}


#' Slop (i.e. extend left/right)
#' @param ranges      data.table(chr, start, end)
#' @param bsgenome    BSgenome, e.g. BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' @param startoffset flank start relative to range start
#' @param endoffset   flank end   relative to range start
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
    startoffset = -22, 
    endoffset   =  22, 
    verbose     = TRUE
){

    # Assert
    assertive.types::assert_is_data.table(ranges)
    assertive.sets::assert_is_subset(c('chr', 'start', 'end'), names(ranges))
    assertive.types::assert_is_a_number(startoffset)
    assertive.types::assert_is_a_number(endoffset)
    assertive.types::assert_is_a_bool(verbose)
    
    # Slop
    start <- end <- NULL
    newranges  <-  data.table::copy(ranges)                 %>% 
                   extract(, start := start + startoffset)  %>%  
                   extract(, end   := end   + endoffset)    %>% 
                   extract()
    newranges [ , chrlength  := GenomeInfoDb::seqlengths(bsgenome)[chr] ]
    tmp <- newranges [, assertive.base::assert_all_are_true(start >= 1)      ]
    tmp <- newranges [, assertive.base::assert_all_are_true(end <= chrlength)]
    newranges[, 'chrlength' := NULL]

    # Return
    if (verbose) cmessage(
                    '\t\t%d slopped ranges: [start%s%d, end%s%d]', 
                    nrow(newranges),
                    csign(startoffset), 
                    abs(startoffset), 
                    csign(endoffset), 
                    abs(endoffset))
    return(newranges)
    
}

reduce <- IRanges::reduce

#' Range merge overlaps
#' @param x   data.table(chr, start, end, strand)
#' @param verbose  logical(1)
#' @return data.table(chr, start, end)
#' @examples 
#' (ranges <- data.table::data.table(
#'             chr    = rep('chr1', 4), 
#'             start  = c(1,    5,   5,  20), 
#'             end    = c(10,  15,  15,  30), 
#'             strand = c('+', '+', '-', '-')))
#' reduce(ranges)
#' @export
setMethod("reduce", signature(x = "data.table"), 
    function(x, verbose = TRUE){
        reduced_ranges <-   dt2gr(x) %>% 
                            GenomicRanges::reduce() %>% 
                            gr2dt()
        if (verbose)  cmessage('\t\t%d ranges after merging overlaps', 
                                nrow(reduced_ranges))
        return(reduced_ranges)
    
    # Alternative data.table solution:    
    # https://stackoverflow.com/a/41748171
    }
)


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
    return(newranges)

}


#' Flank left/right ranges for both strands, merging overlaps
#' @param ranges            data.table(chr, start, end, strand)
#' @param bsgenome          BSgenome, e.g. BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' @param leftstartoffset   numeric(1)
#' @param leftendoffset     numeric(1)
#' @param rightstartoffset  numeric(1)
#' @param rightendoffset    numeric(1)
#' @param verbose           logical(1): report?
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
    leftstartoffset  = -200, 
    leftendoffset    =   -1, 
    rightstartoffset =    1, 
    rightendoffset   =  200, 
    verbose          = TRUE
){
    # Comply
    . <- NULL
    
    # Flank
    if (verbose) cmessage('\tFlank fourways')
    newranges <- rbind(left_flank(
                            ranges,
                            bsgenome    = bsgenome, 
                            startoffset = leftstartoffset, 
                            endoffset   = leftendoffset, 
                            verbose     = verbose
                        ),
                        right_flank(
                            ranges,
                            bsgenome    = bsgenome,
                            startoffset = rightstartoffset,
                            endoffset   = rightendoffset,
                            verbose     = verbose
                        )
                )
    if (verbose) cmessage('\t\t%d ranges combined (left + right)', 
                            nrow(newranges))

    # Complement
    newranges %<>% rbind(complement(., verbose = FALSE))
    newranges %>% data.table::setorderv(c('chr', 'start', 'end'))
    if (verbose) cmessage('\t\t%d ranges after adding strand-complements', 
                            nrow(newranges))

    # Reduce        
    newranges %>% reduce(verbose = verbose) 
}


#' Slop ranges for both strands, merging overlaps
#' @param ranges      data.table(chr, start, end, strand)
#' @param bsgenome    BSgenome, e.g. BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' @param startoffset numeric(1)
#' @param endoffset   numeric(1)
#' @param verbose     logical(1)
#' @return data.table(chr, start, end, strand)
#' @examples
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#' ranges <- read_bed(bedfile)
#' slop_fourways(ranges, bsgebnome)
#' @export
slop_fourways <- function(
    ranges,
    bsgenome,
    startoffset = -22, 
    endoffset   = 22, 
    verbose     = TRUE
){
    # Comply
    . <- NULL
    
    # Slop
    if (verbose) cmessage('\tSlop fourways')
    newranges <-    ranges %>% 
                    slop(
                        bsgenome    = bsgenome,
                        startoffset = startoffset, 
                        endoffset   = endoffset, 
                        verbose     = verbose
                    )
    
    # Complement
    newranges %<>% rbind(complement(., verbose = FALSE))
    if (verbose)   cmessage('\t\t%d ranges after adding strand-complements', 
                            nrow(newranges))
    
    # Reduce
    newranges %>% reduce(verbose = verbose)
}


