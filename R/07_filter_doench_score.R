
#' Add [-4, +3] contextseq
#' 
#' @param cas9s \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome   \code{\link[BSgenome]{BSgenome-class}}
#' @param verbose logical(1)
#' @return character vector
#' @examples 
#' require(magrittr)
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' txdb <- utils::getFromNamespace('TxDb.Mmusculus.UCSC.mm10.knownGene', 
#'                                 'TxDb.Mmusculus.UCSC.mm10.knownGene')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10                                
#' targets  <- bed_to_granges(bedfile, txdb)  %>% 
#'             double_flank() %>% 
#'             add_seq(bsgenome)
#' cas9s    <- find_cas9s(targets)
#' cas9s %<>% add_contextseq(bsgenome)
#' cas9s[1:3]$seq
#' cas9s[1:3]$contextseq
#' @export
add_contextseq <- function(cas9s, bsgenome, verbose = TRUE){
    contexts <- cas9s
    start(contexts) <- start(cas9s) - ifelse(strand(cas9s)=='+', 4, 3)
    end(contexts)   <- end(cas9s)   + ifelse(strand(cas9s)=='+', 3, 4)
    if (verbose) cmessage('\t\tAdd (4-23-3) contextseqs')
    contexts %<>% add_seq(bsgenome, verbose = FALSE)
    cas9s$contextseq <- contexts$seq
    cas9s
}


score_rs1 <- function(
    contextseqs,
    python     = NULL,
    virtualenv = NULL,
    condaenv   = NULL,
    verbose
){
    
    # Assert
    assertive.types::assert_is_character(contextseqs)
    assertive.base::assert_all_are_true(nchar(contextseqs)==30)
    assertive.types::assert_is_a_bool(verbose)
    
    # Message
    if (verbose)  message('\t\tScore contextseqs with ruleset1')
    
    # Score
    CRISPRseek::calculategRNAEfficiency(
        contextseqs,
        baseBeforegRNA       = 4, 
        featureWeightMatrix  =  system.file("extdata", 
                                            "DoenchNBT2014.csv", 
                                            package = "CRISPRseek") %>% 
                                utils::read.csv(header = TRUE)) %>% 
    extract(, 1)
}


score_rs2 <- function(
    contextseqs, 
    python     = NULL, 
    virtualenv = NULL, 
    condaenv   = NULL,
    verbose    = TRUE 
){
    
    # Set python environment
    if (!is.null(python))       reticulate::use_python(python)
    if (!is.null(virtualenv))   reticulate::use_virtualenv(virtualenv)
    if (!is.null(condaenv))     reticulate::use_condaenv(condaenv)
    
    # Assert
    assertive.reflection::assert_is_unix() & 
    assertive.base::is_identical_to_true(
        reticulate::py_module_available('azimuth'))
    assertive.types::assert_is_character(contextseqs)
    assertive.base::assert_all_are_true(nchar(contextseqs)==30)
    assertive.types::assert_is_a_bool(verbose)
    
    # Message
    if (verbose)  message(  '\t\tScore contextseqs with ruleset2',
                            ' (https://github.com/MicrosoftResearch/Azimuth)')
    
    # Score
    score <- contextseq <- NULL
    seqdt   <- data.table::data.table(contextseq = contextseqs)
    scoredt <- data.table::data.table(contextseq = unique(contextseqs))
    azimuth <- reticulate::import("azimuth", delay_load = TRUE)
    numpy   <- reticulate::import("numpy",   delay_load = TRUE)
    scoredt [ , score := azimuth$model_comparison$predict(
                            numpy$array(contextseq)) 
            ]
    seqdt %>% merge(scoredt, by = 'contextseq') %>% extract2('score')
}


#' Filter on Doench score
#' 
#' Score cas9s using Doench ruleset 1 (Doench 2014) or 2 (Doench 2016). 
#' Then filter on set threshold.
#' 
#' ruleset1 is readily available. ruleset2 is accessible after installing 
#' the python module [azimuth](https://github.com/MicrosoftResearch/Azimuth), 
#' and specifying a value for either \code{python} (python binary path), 
#' \code{virtualenv} (python virtual environment dir) or 
#' \code{condaenv} (python conda environment).
#' 
#' @param cas9s     \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome  \code{\link[BSgenome]{BSgenome-class}}
#' @param ruleset    1 (default) or 2 (requires non-NULL value for either python, virtualenv, or condaenv)
#' @param filter     score threshold to filter for
#' @param python     NULL (default) or python binary path with module azimuth
#' @param virtualenv NULL (default) or python virtualenv with module azimuth
#' @param condaenv   NULL (default) or python condaenv with module azimuth
#' @param plot       TRUE (default) or FALSE
#' @param verbose    TRUE (default) or FALSE
#' @return numeric vector
#' @examples
#' # Get cas9s
#'     require(magrittr)
#'     bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'     txdb <- utils::getFromNamespace('TxDb.Mmusculus.UCSC.mm10.knownGene', 
#'                                     'TxDb.Mmusculus.UCSC.mm10.knownGene')
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'     targets <- bed_to_granges(bedfile, txdb)  %>% 
#'                double_flank() %>%
#'                add_seq(bsgenome)
#' 
#' # Find cas9s
#'     cas9s <- find_cas9s(targets)
#'     
#' # Compute and filter on Doench score
#'     cas9s %>% filter_doench_score(bsgenome)
#' 
#' @references 
#' Doench 2014, Rational design of highly active sgRNAs for 
#' CRISPR-Cas9-mediated gene inactivation. Nature Biotechnology,
#' doi: 10.1038/nbt.3026
#' 
#' Doench 2016, Optimized sgRNA design to maximize activity and minimize 
#' off-target effects of CRISPR-Cas9. Nature Biotechnology, 
#' doi: 10.1038/nbt.3437
#' 
#' Python module azimuth: github/MicrosoftResearch/azimuth
#' @export
filter_doench_score <- function(
    cas9s,
    bsgenome, 
    ruleset    = 1,
    filter     = 0.3,
    python     = NULL,
    virtualenv = NULL,
    condaenv   = NULL,
    plot       = TRUE,
    verbose    = TRUE
){
    # Assert
    assertive.types::assert_is_all_of(cas9s, 'GRanges')

    # Add contextseq
    cas9s %<>% add_contextseq(bsgenome, verbose = verbose)
    
    # Score
    scorefun <- switch(ruleset, `1` = score_rs1, `2` = score_rs2)
    cas9s$score <- scorefun(cas9s$contextseq, 
                            verbose     = verbose,
                            python      = python, 
                            virtualenv  = virtualenv, 
                            condaenv    = condaenv)
    
    # Filter
    if (verbose) cmessage('\t\tFilter: score > %s', as.character(filter))
    goodscore <- cas9s %>% extract(cas9s$score > filter)
    
    # Plot
    if (plot){
        grlist <- GRangesList(cas9 = cas9s, goodscore = goodscore)
        names(grlist)[2] <- sprintf('score > %s', as.character(filter))
        plot_karyogram(grlist)
    }
    
    # Message
    if (verbose)  cmessage('\t\tReturn %d/%d cas9seqs at %d/%d cas9 sites: score > %s', 
                           length(unique(goodscore$seq)),
                           length(unique(goodscore$seq)),
                           length(cas9s),
                           length(goodscore), 
                           as.character(filter))
    
    # Return
    cas9s
    
}
