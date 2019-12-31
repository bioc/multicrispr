
#' Add [-4, +3] context
#' 
#' Add context for Doench2016 scoring
#' 
#' @param spacers \code{\link[GenomicRanges]{GRanges-class}}: spacer ranges
#' @param bsgenome   \code{\link[BSgenome]{BSgenome-class}}
#' @param verbose logical(1)
#' @return character vector
#' @examples 
#' # PE example
#' #-----------
#'     require(magrittr)
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#'     gr <- GenomicRanges::GRanges(
#'              seqnames = c(PRNP = 'chr20:4699600',             # snp
#'                           HBB  = 'chr11:5227002',             # snp
#'                           HEXA = 'chr15:72346580-72346583',   # del
#'                           CFTR = 'chr7:117559593-117559595'), # ins
#'              strand   = c(PRNP = '+', HBB = '-', HEXA = '-', CFTR = '+'), 
#'              seqinfo  = BSgenome::seqinfo(bsgenome))
#'     spacers <- find_spacers(extend_for_pe(gr), bsgenome)
#'     add_context(spacers, bsgenome)
#' 
#' # TFBS example
#' #-------------
#'     bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'     targets  <- bed_to_granges(bedfile, 'mm10')  %>% 
#'                 extend() %>% 
#'                 add_seq(bsgenome)
#'     spacers    <- find_spacers(targets, bsgenome)
#'     spacers %<>% add_context(bsgenome)
#'     spacers[1:3]$seq
#'     spacers[1:3]$contextseq
#' @export
add_context <- function(spacers, bsgenome, verbose = TRUE){
    
    # Prevent from stats::start from being used (leads to bug!)
    contexts <- extend(spacers, -4, +6)
    spacers$context <- BSgenome::getSeq(
                            bsgenome, contexts, as.character = TRUE)
    if (verbose) message('\t\tAdd (4-23-3) contextseqs')
    spacers
}


doench2014 <- function(
    contextseqs,
    python     = NULL,
    virtualenv = NULL,
    condaenv   = NULL,
    verbose
){
    
    # Assert
    assert_is_character(contextseqs)
    assert_all_are_true(nchar(contextseqs)==30)
    assert_is_a_bool(verbose)
    
    # Message
    if (verbose)  message('\t\tScore contextseqs with Doench2014')
    
    # Read featureWeightMatrix
    fwfile   <- system.file("extdata/DoenchNBT2014.csv", package = "CRISPRseek")
    fwmatrix <- read.csv(fwfile, header = TRUE)
    
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
    is_identical_to_true(reticulate::py_module_available('azimuth'))
    assert_is_character(contextseqs)
    assert_all_are_true(nchar(contextseqs)==30)
    assert_is_a_bool(verbose)
    
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


#' Score spacers
#' 
#' Score spacers with Doench2014 or Doench2016
#' 
#' Doench2014 is readily available. 
#' Doench2016 is available after installing python module 
#' [azimuth](https://github.com/MicrosoftResearch/Azimuth) (see examples).
#' 
#' @param spacers   \code{\link[GenomicRanges]{GRanges-class}}: spacers
#' @param bsgenome  \code{\link[BSgenome]{BSgenome-class}}
#' @param method   'Doench2014' (default) or 'Doench2016'
#'                 (requires non-NULL argument python, virtualenv, or condaenv)
#' @param python     NULL (default) or python binary path with module azimuth
#' @param virtualenv NULL (default) or python virtualenv with module azimuth
#' @param condaenv   NULL (default) or python condaenv with module azimuth
#' @param verbose    TRUE (default) or FALSE
#' @return numeric vector
#' @examples
#' # PE example
#' #-----------
#'     require(magrittr)
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#'     gr <- GenomicRanges::GRanges(
#'               seqnames = c(PRNP = 'chr20:4699600',             # snp
#'                            HBB  = 'chr11:5227002',             # snp
#'                            HEXA = 'chr15:72346580-72346583',   # del
#'                            CFTR = 'chr7:117559593-117559595'), # ins
#'               strand   = c(PRNP = '+', HBB = '-', HEXA = '-', CFTR = '+'), 
#'               seqinfo  = BSgenome::seqinfo(bsgenome))
#'     spacers <- find_spacers(extend_for_pe(gr), bsgenome)
#'     score_spacers(spacers, bsgenome, 'Doench2014')
#'         # conda create --name azimuthenv python=2.7
#'         # conda activate azimuthenv
#'         # pip install azimuth
#'         # pip install scikit-learn==0.17.1
#'     # score_spacers(spacers, bsgenome, 'Doench2016', condaenv = 'azimuthenv')
#' 
#' # TFBS example
#' #-------------
#'     bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'     gr <- extend(bed_to_granges(bedfile, 'mm10'))
#'     spacers <- find_spacers(gr, bsgenome)
#'     score_spacers(spacers, bsgenome, 'Doench2014')
#'     # score_spacers(spacers, bsgenome, 'Doench2016')
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
score_spacers <- function(
    spacers,
    bsgenome, 
    method     = c('Doench2014', 'Doench2016')[1],
    python     = NULL,
    virtualenv = NULL,
    condaenv   = NULL,
    verbose    = TRUE
){
    # Assert
    assert_is_all_of(spacers, 'GRanges')
    assert_is_a_string(method)
    assert_is_subset(method, c('Doench2014', 'Doench2016'))

    # Add contextseq
    if (verbose)  cmessage('\tScore crispr spacers')
    spacers %<>% add_context(bsgenome, verbose = verbose)
    spacerdt  <- as.data.table(spacers)
    scoredt <- data.table(context = unique(spacerdt$context))
    
    # Score
    scorefun <- switch(method, Doench2014 = doench2014, Doench2016 = doench2016)
    scoredt[ , (method) := scorefun(scoredt$context, verbose=verbose) ]

    # Merge back in and Return
    mergedt  <- merge(spacerdt, scoredt, by='context', sort=FALSE, all.x=TRUE)
    GRanges(mergedt, seqinfo = seqinfo(spacers))
}
