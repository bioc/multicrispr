contextify_start <- function(granges){
    
    # Assert
    assertive.base::assert_is_identical_to_true(methods::is(granges, 'GRanges'))
    
    # Contextify
    GenomicRanges::start(granges) %<>% 
        subtract(ifelse(GenomicRanges::strand(granges)=='+', 4, 3))
    assertive.numbers::assert_all_are_greater_than_or_equal_to(
        GenomicRanges::start(granges), 1)
    
    # Return
    return(granges)
}


contextify_end <- function(granges){

    # Assert
    assertive.base::assert_is_identical_to_true(methods::is(granges, 'GRanges'))
    
    # Contextify
    GenomicRanges::end(granges) %<>% 
        add(ifelse(GenomicRanges::strand(granges) == '+', 3, 4))
    chrlengths  <-  GenomeInfoDb::seqlengths(get_bsgenome(granges)) %>% 
                    extract(as.character(BSgenome::seqnames(granges)))
    assertive.base::assert_all_are_true(granges$contextend < chrlengths)
    
    # Return
    return(granges)
}


score_contextseqs_rs1 <- function(
    contextseqs,
    verbose,
    python     = python,
    virtualenv = virtualenv,
    condaenv   = NULL
){
    
    # Assert
    assertive.types::assert_is_character(contextseqs)
    assertive.base::assert_all_are_true(nchar(contextseqs)==30)
    assertive.types::assert_is_a_bool(verbose)
    
    # Message
    if (verbose)  message('\tScore contextseqs (4-23-3) with ruleset1')
    
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


score_contextseqs_rs2 <- function(
    contextseqs, 
    verbose    = TRUE, 
    python     = NULL, 
    virtualenv = NULL, 
    condaenv   = NULL
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
    if (verbose)  message(  '\tScore contextseqs (4-23-3) with ruleset2',
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


#' @rdname score_cas9ranges
#' @export
score_contextseqs <- function(
    contextseqs, 
    ruleset    = 1, 
    verbose    = TRUE, 
    python     = NULL, 
    virtualenv = NULL, 
    condaenv   = NULL
){
    
    # Assert
    assertive.types::assert_is_a_number(ruleset)
    assertive.sets::assert_is_subset(ruleset, c('1', '2'))
    
    # Score
    scorefun <- switch(ruleset, `1` = score_contextseqs_rs1, 
                                `2` = score_contextseqs_rs2)
    scorefun(
        contextseqs = contextseqs, 
        verbose     = verbose,
        python      = python, 
        virtualenv  = virtualenv, 
        condaenv    = condaenv
    )
    
}

#' Get contextseqs
#' 
#' Get [-3, +3] contextseqs for given cas9ranges
#' 
#' @param cas9ranges GenomicRanges::GRanges
#' @examples 
#' require(magrittr)
#' bedfile <- system.file('extdata/SRF_sites.bed', package = 'multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' cas9ranges  <-  read_bed(bedfile, bsgenome) %>% 
#'                 flank_fourways() %>% 
#'                 find_cas9ranges()
#' cas9ranges[1:3] %>% seqs()
#' cas9ranges[1:3] %>% contextseqs()
#' @export
contextseqs <- function(cas9ranges){
    cas9ranges          %>% 
    contextify_start()  %>%
    contextify_end()    %>% 
    seqs()              %>% 
    as.character()
}


#' Score cas9ranges
#' 
#' Score cas9ranges using ruleset1 (Doench 2014) or ruleset2 (Doench 2016)
#' 
#' ruleset1 is readily available. ruleset2 is accessible after installing 
#' the python module [azimuth](https://github.com/MicrosoftResearch/Azimuth), 
#' and specifying a value for either 'python' (python binary path), 
#' 'virtualenv' (python virtual environment dir) or 
#' 'condaenv' (python conda environment).
#' 
#' @param cas9ranges GenomicRanges::GRanges
#' @param contextseqs character vector with 4-23-3 contextseqs
#' @param ruleset    1 (default) or 2 (only if python module 
#'                   github/MicrosoftResearch/azimuth is installed)
#' @param verbose    logical(1)
#' @param python     NULL (ruleset=1) or path to a python binary (ruleset=2). See details.
#' @param virtualenv NULL (ruleset=1) or directory containing python virtualenv (ruleset=2). See details.
#' @param condaenv   NULL (ruleset=1) or name of condaenv (ruleset=2). See details.
#' @return numeric vector
#' @examples
#' 
#' # Get cas9ranges
#'     require(magrittr)
#'     bedfile <- system.file('extdata/SRF_sites.bed', package = 'multicrispr')
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#'     cas9ranges <- read_bed(bedfile, bsgenome) %>% 
#'                   flank_fourways() %>% 
#'                   find_cas9ranges()
#' # Score
#'     cas9ranges[1:3] %>% score_cas9ranges()
#'     cas9ranges[1:3] %>% contextseqs() %>% score_contextseqs()
#' 
#' @references 
#' Doench 2014, Rational design of highly active sgRNAs for 
#' CRISPR-Cas9-mediated gene inactivation. Nature Biotechnology,
#' doi: 10.1038/nbt.3026
#' 
#' Doench 2016, Optimized sgRNA design to maximize activity and minimize 
#' off-target effects of CRISPR-Cas9. Nature Biotechnology, 
#' doi: 10.1038/nbt.3437
#' @export
score_cas9ranges <- function(
    cas9ranges,
    ruleset    = 1,
    verbose    = TRUE,
    python     = NULL,
    virtualenv = NULL,
    condaenv   = NULL
){
    # Assert
    assertive.base::assert_is_identical_to_true(
        methods::is(cas9ranges, 'GRanges'))

    # Score
    cas9ranges %>% 
    contextseqs() %>%  
    score_contextseqs(
        ruleset    = ruleset, 
        verbose    = verbose, 
        python     = python, 
        virtualenv = virtualenv, 
        condaenv   = condaenv)

}
