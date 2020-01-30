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
    .N <- N <- suffix <- NULL
    dt <- data.table::data.table(x = x)
    dt[, N := .N, by='x']
    dt[N==1, xunique := x]
    dt[N>1, xunique := paste0(x, '_', 1:.N), by = 'x']
    dt[, xunique]
}

#' Make unique names
#' @param x vector
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

