
#' Convert range to sequence
#' @param chr     chromosome vector
#' @param start   start position vector
#' @param end     vector end position vector
#' @param strand  '+|-' vector
#' @param bsgenome  BSgenome object
#' @examples
#' range2seq(
#'     chr      = c('chr2', 'chr1', 'chr17'), 
#'     start    = c(83972889, 85603908, 70963622),
#'     end      = c(83972905, 85603924, 70963638),
#'     strand   = c('+', '-', '-'), 
#'     bsgenome = BSgenome.Mmusculus.UCSC.mm10::Mmusculus)
#' @export 
range2seq <- function(chr, start, end, strand, bsgenome){
    . <- NULL 
    GenomicRanges::GRanges( 
        seqnames = chr,
        ranges   = IRanges::IRanges(start = start, end = end),
        strand   = strand,
        seqinfo  = GenomeInfoDb::seqinfo(bsgenome)) %>%
    BSgenome::getSeq(bsgenome, .) %>%
    as.character()
}


#' Read bedfile as data.table
#' @param bedfile  file path
#' @param verbose  logical(1)
#' @examples
#' bedfile <- system.file('extdata/SRF_sites.bed', package = 'cas9tbflanks')
#' read_bed(bedfile)
#' @importFrom  data.table  :=
#' @export
read_bed <- function(bedfile, verbose = TRUE){

    # Comply
    chr <- start <- end <- strand <- NULL

    # Read
    dt <- data.table::fread(
            bedfile, 
	    select    = c(1:3, 6), 
	    col.names = c('chr', 'start', 'end', 'strand'))

    # Rm strand, since cas9 is strand-agnostic
    dt[, strand := NULL]
    dt %>% data.table::setorder(chr, start, end)

    # Report statistics
    if (verbose){
        dt [ , width := end-start+1]
        dt [ , gap := c(start[2:.N]-end[1:(.N-1)], Inf), by = chr ]
        dt [ , message(sprintf('\t\t%d ranges on %d chromosomes', .N, length(unique(chr)))) ]
        dt [ , message(sprintf('\t\trange width: %s', if (length(unique(width))==1) width[1] else paste0(min(width), ' - ', max(width))))  ]
        dt [ , message(sprintf('\t\tgap   width: %d - %d nt',   min(gap), max(gap[is.finite(gap)]))) ]
        dt [ , c('gap', 'width') := NULL]
    }

    # Return
    return(dt)
}





