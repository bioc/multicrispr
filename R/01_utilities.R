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


#' Get canonical seqlevels
#' @param gr \code{\link[GenomicRanges]{GRanges-class}}
#' @return character vector
#' @seealso \code{\link[GenomeInfoDb]{seqinfo}}
#' @examples
#' 
#' # TxDb mouse
#'     txdbname <- 'TxDb.Mmusculus.UCSC.mm10.knownGene'
#'     txdb <- utils::getFromNamespace(txdbname, txdbname)
#'     gr <- GenomicFeatures::genes(txdb)
#'     GenomeInfoDb::seqlevels(gr)
#'     canonicalseqlevels(gr)
#' 
#' # TxDb Human
#'     txdbname <- 'TxDb.Hsapiens.UCSC.hg38.knownGene'
#'     txdb <- utils::getFromNamespace(txdbname, txdbname)
#'     gr <- GenomicFeatures::genes(txdb)
#'     GenomeInfoDb::seqlevels(gr)
#'     canonicalseqlevels(gr)
#'     
#' # EnsDb mouse
#'     ensdb <- EnsDb.Mmusculus.v98()
#'     gr <- ensembldb::genes(ensdb)
#'     GenomeInfoDb::seqlevels(gr)
#'     canonicalseqlevels(gr)
#'     
#' # EnsDb human
#'     ensdb <- EnsDb.Hsapiens.v98()
#'     gr <- ensembldb::genes(ensdb)
#'     GenomeInfoDb::seqlevels(gr)
#'     canonicalseqlevels(gr)
#'     
#' # Bedfile mouse
#'     bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'     gr <- read_bed(bedfile, 'mm10', plot = FALSE)
#'     GenomeInfoDb::seqlevels(gr)
#'     canonicalseqlevels(gr)
#'     
#' @export
canonicalseqlevels <- function(gr){

    # UCSC
    if ('UCSC' %in% GenomeInfoDb::seqlevelsStyle(gr)){
        GenomeInfoDb::seqlevels(gr) %>% 
        extract(stringi::stri_detect_regex(., '^chr[0-9XYM]+$'))
        
    # NCBI/Ensembl
    } else {
        GenomeInfoDb::seqlevels(gr) %>% 
        extract(stringi::stri_detect_regex(., '^[0-9XYMT]+$'))
    }
    
    #extract(stringi::stri_detect_fixed(., 'chrUn',   negate = TRUE))  %>%
    #extract(stringi::stri_detect_fixed(., '_random', negate = TRUE))  %>%
    #extract(stringi::stri_detect_fixed(., '_fix',    negate = TRUE))  %>%
    #extract(stringi::stri_detect_fixed(., '_alt',    negate = TRUE))
}
