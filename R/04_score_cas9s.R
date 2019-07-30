#' Add cas9 ontarget scores
#' @param cas9dt data.table(chr, cas9start, cas9end, strand)
#' @param bsgenome BSGenome object, e.g. BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' @return data.table(chr, cas9start, cas9end, strand)
#' @examples
#' require(magrittr)
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#' cas9dt <- read_bed(bedfile) %>% slop_fourways() %>% find_cas9s(bsgenome)
#' cas9dt %>% add_contextseq()
#' @export
add_ontargetscores <- function(cas9dt, bsgenome){
    message('\tScore cas9 sites (https://github.com/MicrosoftResearch/Azimuth)')
    cas9dt [ !is.na(cas9start) & strand=='+', contextstart := cas9start-4,                                               ]
    cas9dt [ !is.na(cas9start) & strand=='+', contextend   := cas9end  +3,                                               ]
    cas9dt [ !is.na(cas9start) & strand=='-', contextstart := cas9start-3,                                               ]
    cas9dt [ !is.na(cas9start) & strand=='-', contextend   := cas9end  +4,                                               ]
    cas9dt [ !is.na(cas9start)              , contextseq   := range2seq(chr, contextstart, contextend, strand, bsgenome) ]
    cas9dt [                                , c('contextstart', 'contextend') := NULL                                    ]
    scoredt <- cas9dt [!is.na(cas9start) , .(contextseq = unique(contextseq))  ]
    scoredt [, ontargetscore := azimuth$model_comparison$predict(numpy$array(scoredt$contextseq)) ]
    cas9dt %<>% merge(scoredt, by = 'contextseq', all.x = TRUE)
    cas9dt [ , contextseq := NULL ]
    cas9dt[]
}

