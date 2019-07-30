#' Contextify start/stop
#' @param chr   character vector: chromosome values
#' @param start numeric vector: start positions
#' @param end   numeric vector: end positions
#' @param strand character vector: '+' or '-' values
#' @return numeric vector
#' @examples
#' ranges <- data.table::data.table(
#'             chr    = c('chr1', 'chr1'), 
#'             start  = c(200, 300), 
#'             end    = c(220, 330), 
#'             strand = c('-', '+'))
#' ranges [ , contextstart := contextify_start(chr, start, end, strand) ][]
#' ranges [ , contextend   := contextify_end(  chr, start, end, strand) ][]
#' @export
contextify_start <- function(chr, start, strand){
    tmp <- Reduce( assertive.properties::assert_are_same_length, 
                   list(chr, start, strand))
    assertive.sets::assert_is_subset(unique(strand), c('+', '-'))
    
    ifelse(strand == '+', start - 4, start - 3)
}


#' @rdname contextify_start
#' @export
contextify_end <- function(chr, end, strand){
    tmp <- Reduce( assertive.properties::assert_are_same_length, 
                   list(chr, end, strand))
    assertive.sets::assert_is_subset(unique(strand), c('+', '-'))
    ifelse(strand == '+', end + 3,   end + 4)
}


#' Calculate azimuth ontargetscores
#' @param contextseqs character vector
#' @param verbose logical(1)
#' @return numeric vector
#' @examples
#' if (reticulate::py_module_available('azimuth')){
#'     contextseqs <- c('TGCCCTTATATTGTCTCCAGCAGAAGGTGT',
#'                      'TGCCCTTATATTGTCTCCAGCAGAAGGTGT',
#'                      'CCAAATATTGTCAAGTTGACAACCAGGAAT')
#'     calc_ontargetscore(contextseqs)
#' }
#' @export
score_contextseqs <- function(contextseqs, verbose = TRUE){
    seqdt   <- data.table::data.table(contextseq = contextseqs)
    scoredt <- data.table::data.table(contextseq = unique(contextseqs))
    if (verbose)  message(  '\tScore cas9 sites ', 
                            '(https://github.com/MicrosoftResearch/Azimuth)')
    azimuth <- reticulate::import("azimuth", delay_load = TRUE)
    numpy   <- reticulate::import("numpy",   delay_load = TRUE)
    scoredt [ , score := azimuth$model_comparison$predict(
                            numpy$array(contextseq)) 
            ]
    seqdt %>% merge(scoredt, by = 'contextseq') %>% extract2('score')
}


#' Calc cas9 ontargetscores
#' @param cas9dt data.table(chr, cas9start, cas9end, strand)
#' @param bsgenome BSGenome object, e.g. BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' @return data.table(chr, cas9start, cas9end, strand)
#' @examples
#' if (reticulate::py_module_available('azimuth')){
#'     require(magrittr)
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#'     bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#'     cas9dt <-   read_bed(bedfile)    %>% 
#'                 head(3)              %>% 
#'                 slop_fourways()      %>% 
#'                 find_cas9s(bsgenome)
#'     cas9dt [ , score_cas9ranges(chr, cas9start, cas9end, strand, bsgenome) ]
#'     cas9dt %>% add_ontargetscores(chr, cas9start, cas9end, strand, bsgenome)
#' }
#' @export
score_cas9ranges <- function(
    chr, 
    start, 
    end, 
    strand, 
    bsgenome, 
    verbose = TRUE
){
    range2seq(  chr, 
                contextify_start(chr, start, strand), 
                contextify_end(  chr, end,   strand), 
                strand, 
                bsgenome) %>% 
    score_contextseqs(verbose)

}

