cduplicated <- function(x)  duplicated(x) | duplicated(x, fromLast = TRUE)

cmessage <- function(...) message(sprintf(...))


num2scalarstr <- function(x){
    if (length(unique(x))==1){
        x[1]
    } else {
        x %<>% magrittr::extract(!is.na(x))
        x %<>% magrittr::extract(is.finite(x))
        paste0(min(x), ' - ', max(x))
    }
}

csign <- function(x) if (sign(x)==-1) '-' else '+'

#' Switch between data.table and GenomicRanges::GRanges representation
#' @param dt data.table
#' @param gr GenomicRanges::GRanges
#' @return IRanges
#' @examples
#' require(magrittr) 
#' dt <- data.table::data.table(
#'     chr    = 'chr1',
#'     start  = c( 1,  5, 5), 
#'     end    = c(10, 15, 15), 
#'     strand = c('+', '+', '-')
#' )
#' dt %>% dt2gr()
#' dt %>% dt2gr() %>% gr2dt()
#' @export
dt2gr <- function(dt){
    
    # Assert
    assertive.types::is_data.table(dt)
    assertive.sets::assert_is_subset(
        c('chr', 'start', 'end', 'strand'), names(dt))
    
    # Convert
    GenomicRanges::GRanges(
        seqnames = dt$chr, 
        ranges   = IRanges::IRanges(start = dt$start, end = dt$end), 
        strand   = dt$strand
    )
}


#' @rdname dt2gr
#' @export
gr2dt <- function(gr){
    
    # Assert
    assertive.base::assert_is_identical_to_true(methods::is(gr, 'GRanges'))
    
    # Convert
    data.table::data.table(
        chr    = gr %>% GenomicRanges::seqnames()  %>%  as.vector(), 
        start  = gr %>% GenomicRanges::start(),
        end    = gr %>% GenomicRanges::end(),
        strand = gr %>% GenomicRanges::strand()    %>% as.vector()
    )
}