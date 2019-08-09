cduplicated <- function(x)  duplicated(x) | duplicated(x, fromLast = TRUE)

#' Message with sprintf syntax
#' @param ... passed to message(sprintf(...))
#' @return NULL
#' @export
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

