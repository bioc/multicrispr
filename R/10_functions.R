



#' Convert range data.table into sequence data.table
#' @param rangedt data.table(chr, start, end)
#' @param verbose logical(1)
#' @return data.table
#' @examples 
#' rangedt <- data.table::data.table(
#'     chr      = c('chr2', 'chr1', 'chr17'),
#'     start    = c(83972889, 85603908, 70963622),
#'     end      = c(83972905, 85603924, 70963638),
#'     strand   = c('+', '-', '-'))
#'     
#' rangedt2seqdt(rangedt, BSgenome.Mmusculus.UCSC.mm10::Mmusculus) 
rangedt2seqdt <- function(rangedt, bsgenome, verbose = TRUE){
    minseqs  <- rangedt [ , range2seq(chr, start, end, strand = '-', bsgenome)]
    plusseqs <- rangedt [ , range2seq(chr, start, end, strand = '+', bsgenome)]
    seqdt <- data.table::data.table(seq = c(minseqs, plusseqs))
    if (verbose) cmessage('\t\t%d strand-specific flank seqs', nrow(seqdt))
    seqdt %<>% extract(, .(nflanks = .N), by = seq)
    if (verbose){
        cmessage(
            '\t\tRetain %d unique seqs (after removing %d replicated seqs)',
            nrow(seqdt),
            seqdt[ nflanks>1, sum(nflanks) - .N])
    }
    return(seqdt)
}
