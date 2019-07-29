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
#' @param startoffset flank start relative to range start
#' @param endoffset   flank end   relative to range start
#' @param verbose     logical(1)
#' @return data.table
#' @export
#' @examples 
#' require(magrittr)
#' bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#' ranges <- read_bed(bedfile)
#' ranges %>% head(1)
#' ranges %>% head(1) %>% left_flank()
left_flank <- function(
    ranges, 
    startoffset = -200, 
    endoffset   = -1, 
    verbose     = TRUE
){
    start <- end <- NULL
    flankranges <-  data.table::copy(ranges) %>% 
                    extract(, end   := start + endoffset)   %>%  # dont switch
                    extract(, start := start + startoffset) %>%  # lines!
                    extract()
    
    cmessage('\t\t%d left  flanks : [start%s%d, start%s%d]', 
            nrow(flankranges),
            csign(startoffset), 
            abs(startoffset), 
            csign(endoffset),
            abs(endoffset))

    return(flankranges)
}


#' Right flank 
#' @param ranges      data.table(chr, start, end)
#' @param startoffset flank start relative to range start
#' @param endoffset   flank end   relative to range start
#' @param verbose     logical(1)
#' @return data.table
#' @export
#' @examples 
#' require(magrittr)
#' bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#' ranges <- read_bed(bedfile)
#' ranges %>% head(1)
#' ranges %>% head(1) %>% right_flank()
#' @export
right_flank <- function(
    ranges, 
    startoffset = 1, 
    endoffset   = 200, 
    verbose     = TRUE
){
    start <- end <- NULL
    flankranges <-  data.table::copy(ranges)               %>% 
                    extract(, start := end + startoffset)  %>% 
                    extract(, end   := end + endoffset)    %>% 
                    extract()
    
    cmessage('\t\t%d right flanks : [end%s%d, end%s%d]', 
            nrow(flankranges),
            csign(startoffset), 
            abs(startoffset), 
            csign(endoffset), 
            abs(endoffset))

    return(flankranges)
}


#' Slop (i.e. extend left/right)
#' @param ranges      data.table(chr, start, end)
#' @param startoffset flank start relative to range start
#' @param endoffset   flank end   relative to range start
#' @param verbose     logical(1)
#' @return data.table
#' @export
#' @examples 
#' require(magrittr)
#' bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#' ranges <- read_bed(bedfile)
#' ranges %>% head(1)
#' ranges %>% head(1) %>% slop()
#' @export
slop <- function(
    ranges, 
    startoffset = -22, 
    endoffset   =  22, 
    verbose     = TRUE
){

    start <- end <- NULL
    
    expandedranges <-   data.table::copy(ranges)                 %>% 
                        extract(, start := start + startoffset)  %>%  
                        extract(, end   := end   + endoffset)    %>% 
                        extract()

    if (verbose) cmessage(
                    '\t\t%d slopped ranges: [start%s%d, end%s%d]', 
                    nrow(expandedranges),
                    csign(startoffset), 
                    abs(startoffset), 
                    csign(endoffset), 
                    abs(endoffset))
    
    return(expandedranges)
    
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
    strand <- NULL
    assertive.sets::assert_is_subset(unique(ranges$strand), c('-', '+'))

    complstrand <-  data.table::copy(ranges) %>% 
                    extract(, strand := ifelse(strand=='+', '-', '+')) %>% 
                    extract()
    
    if (verbose) cmessage('\t\t%d strand-complementary ranges', 
                            nrow(complstrand))
    return(complstrand)

}


#' Flank left/right ranges for both strands, merging overlaps
#' 
#' @param ranges            data.table(chr, start, end, strand)
#' @param leftstartoffset   numeric(1)
#' @param leftendoffset     numeric(1)
#' @param rightstartoffset  numeric(1)
#' @param rightendoffset    numeric(1)
#' @param verbose           logical(1): report?
#' @return data.table(chr, start, end, strand)
#' @examples 
#' bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#' ranges <- read_bed(bedfile)
#' flank_fourways(ranges)
##' @export
flank_fourways <- function(
    ranges,
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
    flankranges <- rbind(left_flank(ranges, 
                                    startoffset = leftstartoffset, 
                                    endoffset   = leftendoffset, 
                                    verbose     = verbose),
                        right_flank(ranges,
                                    startoffset = rightstartoffset,
                                    endoffset   = rightendoffset,
                                    verbose = verbose) )
    if (verbose) cmessage('\t\t%d ranges combined (left + right)', 
                            nrow(flankranges))

    # Complement
    flankranges %<>% rbind(complement(., verbose = FALSE))
    flankranges %>% data.table::setorderv(c('chr', 'start', 'end'))
    if (verbose) cmessage('\t\t%d ranges after adding strand-complements', 
                            nrow(flankranges))

    # Reduce        
    flankranges %>% reduce(verbose = verbose) 
}


#' Slop ranges for both strands, merging overlaps
#' @param ranges data.table(chr, start, end, strand)
#' @param startoffset numeric(1)
#' @param endoffset   numeric(1)
#' @param verbose     logical(1)
#' @return data.table(chr, start, end, strand)
#' @examples
#' bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#' ranges <- read_bed(bedfile)
#' slop_fourways(ranges)
#' @export
slop_fourways <- function(
    ranges,
    startoffset = -22, 
    endoffset   = 22, 
    verbose     = TRUE
){
    # Comply
    . <- NULL
    
    # Slop
    if (verbose)   cmessage('\tSlop fourways')
    sloppedranges  <-   ranges %>% 
                        slop(
                            startoffset = startoffset, 
                            endoffset   = endoffset, 
                            verbose     = verbose
                        )
    
    # Complement
    sloppedranges %<>% rbind(complement(., verbose = FALSE))
    if (verbose)   cmessage('\t\t%d ranges after adding strand-complements', 
                            nrow(sloppedranges))
    
    # Reduce
    sloppedranges %>% reduce(verbose = verbose)
}


