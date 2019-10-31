
#' Add [-4, +3] contextseq
#' 
#' @param cas9s \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome   \code{\link[BSgenome]{BSgenome-class}}
#' @param verbose logical(1)
#' @return character vector
#' @examples 
#' require(magrittr)
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' targets  <- bed_to_granges(bedfile, 'mm10')  %>% 
#'             double_flank() %>% 
#'             add_seq(bsgenome)
#' cas9s    <- find_cas9s(targets)
#' cas9s %<>% add_contextseq(bsgenome)
#' cas9s[1:3]$seq
#' cas9s[1:3]$contextseq
#' @export
add_contextseq <- function(cas9s, bsgenome, verbose = TRUE){
    
    # Prevent from stats::start from being used (leads to bug!)
    start  <- GenomicRanges::start;     `start<-`  <- GenomicRanges::`start<-`
    end    <- GenomicRanges::end;       `end<-`    <- GenomicRanges::`end<-`
    strand <- GenomicRanges::strand
    
    contexts <- cas9s
    start(contexts) <- start(cas9s) - ifelse(strand(cas9s)=='+', 4, 3)
    end(contexts)   <- end(cas9s)   + ifelse(strand(cas9s)=='+', 3, 4)
    if (verbose) cmessage('\t\tAdd (4-23-3) contextseqs')
    contexts %<>% add_seq(bsgenome, verbose = FALSE)
    cas9s$contextseq <- contexts$seq
    cas9s
}


doench2014 <- function(
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
    if (verbose)  message('\t\tScore contextseqs with Doench2014')
    
    # Read featureWeightMatrix
    fwfile   <- system.file("extdata/DoenchNBT2014.csv", package = "CRISPRseek")
    fwmatrix <- utils::read.csv(fwfile, header = TRUE)
    
    # Score and return
    CRISPRseek::calculategRNAEfficiency(contextseqs,
                                        baseBeforegRNA       = 4, 
                                        featureWeightMatrix  =  fwmatrix) %>% 
    magrittr::extract(, 1)
}


doench2016 <- function(
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
    assertive.base::is_identical_to_true(
        reticulate::py_module_available('azimuth'))
    assertive.types::assert_is_character(contextseqs)
    assertive.base::assert_all_are_true(nchar(contextseqs)==30)
    assertive.types::assert_is_a_bool(verbose)
    
    # Message
    if (verbose)  message(  '\t\tScore contextseqs with Doench2016',
                            ' (https://github.com/MicrosoftResearch/Azimuth)')
    
    # Score
    azi <- reticulate::import('azimuth.model_comparison', delay_load = TRUE)
    azi$predict(
        reticulate::np_array(contextseqs), 
        aa_cut                 = NULL, 
        percent_peptide        = NULL, 
        model                  = NULL, 
        model_file             = NULL, 
        pam_audit              = TRUE, 
        length_audit           = TRUE, 
        learn_options_override = NULL)
}


#' Score cas9 sites
#' 
#' Score cas9s with Doench2014 or Doench2016 model.
#' 
#' Doench2014 is readily available. 
#' Doench2016 is available after installing python module 
#' [azimuth](https://github.com/MicrosoftResearch/Azimuth) (see examples).
#' 
#' @param cas9s     \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome  \code{\link[BSgenome]{BSgenome-class}}
#' @param method   'Doench2014' (default) or 'Doench2016'
#'                 (requires non-NULL argument python, virtualenv, or condaenv)
#' @param python     NULL (default) or python binary path with module azimuth
#' @param virtualenv NULL (default) or python virtualenv with module azimuth
#' @param condaenv   NULL (default) or python condaenv with module azimuth
#' @param verbose    TRUE (default) or FALSE
#' @return numeric vector
#' @examples
#' # Define targets
#'     require(magrittr)
#'     bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'     targets <- bed_to_granges(bedfile, 'mm10')  %>% 
#'                double_flank() %>%
#'                add_seq(bsgenome)
#' 
#' # Find cas9s
#'     cas9s <- find_cas9s(targets)
#'     
#' # Score with Doench2014
#'     score_cas9s(cas9s[1:10], bsgenome)
#'         
#' # Score with Doench2016
#'     # First install python module azimuth, perhaps in a conda env:
#'         # install conda for python 2.7
#'         # conda create --name azimuthenv python=2.7
#'         # conda activate azimuthenv
#'         # pip install azimuth
#'         # pip install scikit-learn==0.17.1
#'         
#'     # Then call score_cas9s with reference to conda env
#'         # score_cas9s(cas9s[1:10], bsgenome, 'Doench2016', 
#'         #             condaenv = 'azimuthenv')
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
score_cas9s <- function(
    cas9s,
    bsgenome, 
    method     = c('Doench2014','Doench2016')[1],
    python     = NULL,
    virtualenv = NULL,
    condaenv   = NULL,
    verbose    = TRUE
){
    # Assert
    assertive.types::assert_is_all_of(cas9s, 'GRanges')

    # Add contextseq
    if (verbose)  cmessage('\tScore cas9s')
    cas9s %<>% add_contextseq(bsgenome, verbose = verbose)
    cas9dt  <- data.table::as.data.table(cas9s)
    scoredt <- data.table::data.table(contextseq = unique(cas9dt$contextseq))
    
    # Score
    scorefun <- switch(method, Doench2014 = doench2014, Doench2016 = doench2016)
    scoredt[ , (method) := scorefun(scoredt$contextseq, verbose=verbose) ]

    # Merge back in and Return
    cas9smerged  <- merge(  cas9dt, scoredt, by = 'contextseq', sort = FALSE, 
                            all.x = TRUE) %>%
                    methods::as('GRanges')
    GenomeInfoDb::seqinfo(cas9smerged) <- GenomeInfoDb::seqinfo(cas9s)
    cas9smerged
}
