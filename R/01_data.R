#' Get EnsDb.Mmusculus.v98 from AnnotationHub
#' @return \code{\link[ensembldb]{EnsDb-class}}
#' @export
#' @examples 
#' EnsDb.Mmusculus.v98()
EnsDb.Mmusculus.v98 <- function(){
    hub <- AnnotationHub::AnnotationHub()
    #hub <- hub[hub$species == 'Mus musculus' & hub$rdataclass == 'EnsDb']
    #sort(hub$title)
    #AnnotationHub::query(hub, '98') # 'AH75036'
    hub[["AH75036"]]
}

#' Get EnsDb.Hsapiens.v98 from AnnotationHub
#' @return \code{\link[ensembldb]{EnsDb-class}}
#' @export
#' @examples 
#' EnsDb.Hsapiens.v99()
EnsDb.Hsapiens.v99 <- function(){
    hub <- AnnotationHub::AnnotationHub()
    #hub <- hub[hub$species == 'Homo sapiens' & hub$rdataclass == 'EnsDb']
    #sort(hub$title)
    #AnnotationHub::query(hub, '98') # 'AH75011'
    #AnnotationHub::query(hub, '99') # 'AH78783'
    hub[["AH78783"]]
}
