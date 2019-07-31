#' Contextify start/stop
#' @param chr   character vector: chromosome values
#' @param start numeric vector: start positions
#' @param end   numeric vector: end positions
#' @param strand character vector: '+' or '-' values
#' @param bsgenome BSGenome, e.g. BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' @return numeric vector
#' @examples
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' contextify_start('chr1', 200,  '+', bsgenome)
#' contextify_end(  'chr1', 200,  '-', bsgenome)
#' @export
contextify_start <- function(chr, start, strand, bsgenome){
    
    # Assert
    assertive.types::assert_is_character(chr)
    assertive.types::assert_is_numeric(start)
    assertive.types::assert_is_character(strand)
    assertive.sets::assert_is_subset(unique(strand), c('+', '-'))
    tmp <- Reduce( assertive.properties::assert_are_same_length, 
                   list(chr, start, strand))
    assertive.base::assert_is_identical_to_true(
        methods::is(bsgenome, 'BSgenome'))
    
    # Contextify
    contextstart <- ifelse(strand == '+', start - 4, start - 3)
    contextstart %>% assertive.numbers::assert_all_are_in_closed_range(
                        1, 
                        GenomeInfoDb::seqlengths(bsgenome)[[chr]])
    
    # Return
    return(contextstart)
}


#' @rdname contextify_start
#' @export
contextify_end <- function(chr, end, strand, bsgenome){

    # Assert
    assertive.types::assert_is_character(chr)
    assertive.types::assert_is_numeric(end)
    assertive.types::assert_is_character(strand)
    assertive.sets::assert_is_subset(unique(strand), c('+', '-'))
    tmp <- Reduce( assertive.properties::assert_are_same_length, 
                   list(chr, end, strand))
    assertive.base::assert_is_identical_to_true(
        methods::is(bsgenome, 'BSgenome'))
    
    # Contextify
    contextend <- ifelse(strand == '+', end + 3,   end + 4)
    contextend %>% assertive.numbers::assert_all_are_in_closed_range(
                        1, 
                        GenomeInfoDb::seqlengths(bsgenome)[[chr]])
    
    # Return
    return(contextend)
}


#' Score contextseqs
#' @param contextseqs character vector
#' @param verbose logical(1)
#' @return numeric vector
#' @examples
#' \dontrun{
#' if (reticulate::py_module_available('azimuth')){
#'     contextseqs <- c('TGCCCTTATATTGTCTCCAGCAGAAGGTGT',
#'                      'TGCCCTTATATTGTCTCCAGCAGAAGGTGT',
#'                      'CCAAATATTGTCAAGTTGACAACCAGGAAT')
#'     score_contextseqs(contextseqs)
#' }
#' }
#' @export
score_contextseqs <- function(contextseqs, verbose = TRUE){
    
    # Assert
    tmp <- Reduce(assertive.base::assert_are_identical, nchar(contextseqs))
    assertive.types::assert_is_a_bool(verbose)
    
    # Score
    score <- contextseq <- NULL
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


#' Score cas9ranges
#' @param chr       character vector: chromosome values
#' @param start     numeric vector: start positions
#' @param end       numeric vector: end positions
#' @param strand    character vector: '+' or '-' values
#' @param bsgenome  BSGenome, e.g. BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' @param verbose   logical(1)
#' @return numeric vector
#' @examples
#' \dontrun{
#' if (reticulate::py_module_available('azimuth')){
#'     require(magrittr)
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#'     bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#'     cas9dt <-   read_bed(bedfile)    %>% 
#'                 head(3)              %>% 
#'                 slop_fourways()      %>% 
#'                 find_cas9s(bsgenome)
#'     cas9dt [ , score_cas9ranges(chr, cas9start, cas9end, strand, bsgenome) ]
#'     cas9dt %>% score_cas9ranges(chr, cas9start, cas9end, strand, bsgenome)
#' }
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
                contextify_start(chr, start, strand, bsgenome), 
                contextify_end(  chr, end,   strand, bsgenome), 
                strand, 
                bsgenome) %>% 
    score_contextseqs(verbose)

}

