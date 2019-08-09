
#=============================================================================
# bedfile -> data.table -> GRanges
#=============================================================================



#' Get BSgenome
#' @param granges GenomicRanges::GRanges
#' @return BSgenome
#' @examples 
#' require(magrittr)
#' bedfile  <- system.file('extdata/SRF_sites.bed', package='crisprapex')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' granges <- read_bed(bedfile, bsgenome) %>% flank_fourways()
#' granges %>% get_bsgenome()
#' @export
get_bsgenome <- function(granges){
    
    assertive.base::assert_is_identical_to_true(is(granges, 'GRanges'))
    genome <- GenomeInfoDb::genome(granges) %>% unname() %>% unique()
    assertive.types::assert_is_a_string(genome)
    BSgenome::getBSgenome(genome)
    
}

#' Assert that object is a genomic ranges datatable
#' @param x data.table
#' @param bsgenome BSgenome, e.g. BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' @return invisible(x)
#' @examples 
#' require(magrittr)
#' bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' x <- read_bed(bedfile, bsgenome) %>% data.table::as.data.table()
#' assert_is_granges_datatable(x, bsgenome)
#' @export
assert_is_granges_datatable <- function (x, bsgenome){
    
    x        %>% assertive.types::assert_is_data.table()
    names(x) %>% assertive.sets::assert_is_superset(
                    c('seqnames', 'start', 'end', 'strand'))
    x$chr    %>% assertive.sets::assert_is_subset(BSgenome::seqnames(bsgenome))
    x$strand %>% assertive.sets::assert_is_subset(c('-', '+'))
    x$start  %>% assertive.numbers::assert_all_are_greater_than_or_equal_to(1)
    x$end    %>% assertive.numbers::assert_all_are_less_than_or_equal_to(
                    GenomeInfoDb::seqlengths(bsgenome)[x$chr])
    invisible(x)
}



#' Convert data.table into GenomicRanges::GRanges
#' @param x data.table
#' @param bsgenome BSgenome, e.g. BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' @return GenomicRanges::GRanges
#' @examples
#' require(magrittr)
#' bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' x <- read_bed(bedfile, bsgenome)
#' x %>% data.table::as.data.table() %>% as.granges(bsgenome = bsgenome)
#' @export
as.granges <- function(dt, bsgenome){
    
    # Assert
    assert_is_granges_datatable(dt, bsgenome)

    # Convert
    GenomicRanges::GRanges(
        seqnames = dt$seqnames, 
        ranges   = IRanges::IRanges(start = dt$start, end = dt$end), 
        strand   = dt$strand, 
        seqinfo  = BSgenome::seqinfo(bsgenome)
    )

}



#' Read bedfile as data.table
#' @param bedfile  file path
#' @param bsgenome BSgenome, e.g. BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' @param verbose  logical(1)
#' @param rm_duplicates  logical(1)
#' @return data.table(chr, start, end, strand) 
#' @examples
#' bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' read_bed(bedfile, bsgenome)
#' @importFrom  data.table  :=
#' @export
read_bed <- function(
    bedfile, 
    bsgenome,
    verbose       = TRUE, 
    rm_duplicates = TRUE
){
    # Assert
    assertive.files::assert_all_are_existing_files(bedfile)
    assertive.types::assert_is_a_bool(verbose)
    assertive.types::assert_is_a_bool(rm_duplicates)

    # Comply
    seqnames <- start <- end <- strand <- .N <- gap <- width <- NULL

    # Read
    if (verbose) cmessage('\tRead %s', bedfile)
    dt <- data.table::fread(
            bedfile,
            select    = c(seq_len(3), 6),
            col.names = c('seqnames', 'start', 'end', 'strand'))
    dt %>% data.table::setorderv(c('seqnames', 'start', 'end', 'strand'))
    if (verbose) cmessage('\t%d ranges on %d chromosomes', 
                            nrow(dt), length(unique(dt$seqnames)))
    
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
        dt [ , gap := c(start[2:.N]-end[seq_len(.N-1)], Inf), by = seqnames ]
        cmessage('\t\t%s NT wide', num2scalarstr(dt$width))
        cmessage('\t\t%s NT apart', num2scalarstr(dt$gap))
        dt [ , c('gap', 'width') := NULL]
    }
    
    # Return
    dt %>% as.granges(bsgenome)
     
}

# Get sequence values
# @param granges GenomicRanges::GRanges
# @param as_character logical(1)
# @return DNAStringSet or character vector
# @examples
# bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
# bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
# granges <- read_bed(bedfile, bsgenome)
# sequence(granges)
# sequence(granges, as_character = FALSE)
# @export
# sequence <- function(granges, as_character = FALSE){
#     dnastringset <- BSgenome::getSeq(get_bsgenome(granges), granges)
#     if (as_character) return(as.character(dnastringset))
#     return(dnastringset)
# }


