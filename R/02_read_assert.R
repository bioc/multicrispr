
#=============================================================================
# data.table -> GRanges
#=============================================================================

#' Convert data.table into GenomicRanges::GRanges
#' @param x data.table
#' @param bsgenome BSgenome, e.g. BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' @return GenomicRanges::GRanges
#' @examples
#' require(magrittr)
#' bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#' x <- read_bed(bedfile)
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' x %>% as.granges(bsgenome = bsgenome)
#' @export
as.granges <- function(dt, bsgenome){
    
    # Assert
    assert_is_granges_datatable(dt, bsgenome)

    # Convert
    GenomicRanges::GRanges(
        seqnames = dt$chr, 
        ranges   = IRanges::IRanges(start = dt$start, end = dt$end), 
        strand   = dt$strand, 
        seqinfo  = BSgenome::seqinfo(bsgenome)
    )
}



#=============================================================================
# GRanges -> data.table
#=============================================================================

#' @export
as.data.table <- data.table::as.data.table

#' Convert GenomicRanges into data.table
#' @param x 
#' @return data.table(chr, start, end, end, strand)
#' @examples 
#' require(magrittr)
#' x <- read_bed(system.file('extdata/SRF_sites.bed', package = 'crisprapex'))
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' x %>% as.granges(bsgenome = bsgenome) %>% as.data.table()
#' @export
as.data.table.GRanges <- function(x){
    
    # Assert
    assertive.base::assert_is_identical_to_true(methods::is(x, 'GRanges'))
    
    # Convert
    data.table::data.table(
        chr    = x %>% GenomicRanges::seqnames()  %>%  as.vector(), 
        start  = x %>% GenomicRanges::start(),
        end    = x %>% GenomicRanges::end(),
        strand = x %>% GenomicRanges::strand()    %>% as.vector()
    )
}

#=============================================================================
# Read into GRanges
#=============================================================================


#' Read bedfile as data.table
#' @param bedfile        file path
#' @param verbose        logical(1)
#' @param rm_duplicates  logical(1)
#' @return data.table(chr, start, end, strand) 
#' @examples
#' bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#' read_bed(bedfile)
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
    dt
}


#=============================================================================
# Assert
#=============================================================================


#' Assert that object is a genomic ranges datatable
#' @param x data.table
#' @param bsgenome BSgenome, e.g. BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' @return invisible(x)
#' @examples 
#' x <- read_bed(system.file('extdata/SRF_sites.bed', package = 'crisprapex'))
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' assert_is_granges_datatable(x, bsgenome)
#' @export
assert_is_granges_datatable <- function (x, bsgenome){
    
    x        %>% assertive.types::assert_is_data.table()
    names(x) %>% assertive.sets::assert_is_superset(
                    c('chr', 'start', 'end', 'strand'))
    x$chr    %>% assertive.sets::assert_is_subset(BSgenome::seqnames(bsgenome))
    x$strand %>% assertive.sets::assert_is_subset(c('-', '+'))
    x$start  %>% assertive.numbers::assert_all_are_greater_than_or_equal_to(1)
    x$end    %>% assertive.numbers::assert_all_are_less_than_or_equal_to(
                    GenomeInfoDb::seqlengths(bsgenome)[x$chr])
    invisible(x)
}

