cduplicated <- function(x)  duplicated(x) | duplicated(x, fromLast = TRUE)

cmessage <- function(...) message(sprintf(...))


num2scalarstr <- function(x){
    if (length(unique(x))==1){
        x[1]
    } else {
        x %<>% extract(!is.na(x))
        x %<>% extract(is.finite(x))
        paste0(min(x), ' - ', max(x))
    }
}

csign <- function(x) if (sign(x)==-1) '-' else '+'

uniquify <- function(x){
    .N <- N <- suffix <- xunique <- NULL
    dt <- data.table::data.table(x = x)
    dt[, N := .N, by='x']
    dt[N==1, xunique := x]
    dt[N>1, xunique := paste0(x, '_', seq_len(.N)), by = 'x']
    dt[, xunique]
}

#' Make unique names
#' @param x vector
#' @param prefix string: prefix with which to start names
#' @return character vector with unique names
make_unique_names <- function(x, prefix='T'){
    
    if (has_names(x)) return(uniquify(names(x)))
    
    paste0(prefix, formatC(seq_along(x), 
                                digits = floor(log10(length(x))), 
                                flag = 0))
}

name_uniquely <- function(gr, prefix = 'x'){
    names(gr) <- make_unique_names(gr, prefix)
    gr
}

#' GRanges <-> data.table
#' 
#' @param gr      \code{\link[GenomicRanges]{GRanges-class}}
#' @param dt      data.table
#' @param seqinfo \code{\link[GenomeInfoDb]{Seqinfo-class}}
#' @aliases dt2gr
#' @return data.table (gr2dt) or GRanges (dt2gr)
#' @examples
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#' gr <- char_to_granges(c(PRNP = 'chr20:4699600:+',             # snp
#'                             HBB  = 'chr11:5227002:-',             # snp
#'                             HEXA = 'chr15:72346580-72346583:-',   # del
#'                             CFTR = 'chr7:117559593-117559595:+'), # ins
#'                           bsgenome)
#' (dt <- gr2dt(gr))
#' (gr <- dt2gr(dt, BSgenome::seqinfo(bsgenome)))
#' @export
gr2dt <- function(gr){
    dt <- as.data.table(gr)
    if (has_names(gr)){
        dt[ , names := names(gr)]
    }
    dt[]
}

#' @rdname gr2dt
#' @export
dt2gr <- function(dt, seqinfo){
    gr <- GRanges(dt, seqinfo = seqinfo)
    if ('names' %in% names(mcols(gr))){
        names(gr) <- gr$names
        gr$names <- NULL
    }
    gr
}

