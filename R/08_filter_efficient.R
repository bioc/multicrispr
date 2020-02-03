
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
#'     gr <- char_to_granges(c(PRNP = 'chr20:4699600:+',             # snp
#'                             HBB  = 'chr11:5227002:-',             # snp
#'                             HEXA = 'chr15:72346580-72346583:-',   # del
#'                             CFTR = 'chr7:117559593-117559595:+'), # ins
#'                           bsgenome)
#'     spacers <- find_spacers(extend_for_pe(gr), bsgenome, complement=FALSE)
#'     (add_context(spacers, bsgenome))
#' 
#' # TFBS example
#' #-------------
#'     bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'     targets  <- extend(bed_to_granges(bedfile, 'mm10'))
#'     spacers    <- find_spacers(targets, bsgenome)
#'     (spacers %<>% add_context(bsgenome))
#' @export
add_context <- function(spacers, bsgenome, verbose = TRUE){
    
    # Prevent from stats::start from being used (leads to bug!)
    contexts <- extend(spacers, -4, +6)
    spacers$crisprcontext <- BSgenome::getSeq(
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


#' Add efficiency scores
#' 
#' Add Doench2014 or Doench2016 efficiency scores
#' 
#' @param spacers  \code{\link[GenomicRanges]{GRanges-class}}: spacers
#' @param bsgenome \code{\link[BSgenome]{BSgenome-class}}
#' @param method   'Doench2014' (default) or 'Doench2016'
#'                 (requires non-NULL argument python, virtualenv, or condaenv)
#' @param python     NULL (default) or python binary path with module azimuth
#' @param virtualenv NULL (default) or python virtualenv with module azimuth
#' @param condaenv   NULL (default) or python condaenv with module azimuth
#' @param verbose    TRUE (default) or FALSE
#' @param plot       TRUE (default) or FALSE
#' @param alpha_var  NULL or string: var mapped to alpha in plot
#' @return numeric vector
#' @examples
#' # PE example
#' #-----------
#'     require(magrittr)
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#'     targets <- char_to_granges(c(PRNP = 'chr20:4699600:+',             # snp
#'                                  HBB  = 'chr11:5227002:-',             # snp
#'                                  HEXA = 'chr15:72346580-72346583:-',   # del
#'                                  CFTR = 'chr7:117559593-117559595:+'), # ins
#'                                bsgenome)
#'     spacers <- find_pe_spacers(targets, bsgenome)
#'    #spacers <- find_spacers(extend_for_pe(gr), bsgenome, complement = FALSE)
#'     (spacers %<>% add_efficiency(bsgenome, 'Doench2014'))
#'         # conda create --name azimuthenv python=2.7
#'         # conda activate azimuthenv
#'         # pip install azimuth
#'         # pip install scikit-learn==0.17.1
#'     # spacers %<>% add_efficiency(bsgenome, 'Doench2016', condaenv = 'azimuthenv')
#'     # filter_efficient(spacers, bsgenome, 'Doench2016', 0.4, condaenv='azimuthenv')
#' # TFBS example
#' #-------------
#'     bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'     targets <- extend(bed_to_granges(bedfile, 'mm10'))
#'     spacers <- find_spacers(targets, bsgenome)
#'     # (spacers %<>% add_specificity(targets, bsgenome))
#'     # (spacers %>% add_efficiency(bsgenome, 'Doench2014'))
#'     # (spacers %>% add_efficiency(bsgenome, 'Doench2016'))
#'     # spacers %>% filter_efficient(bsgenome, 'Doench2016', 0.4)
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
add_efficiency <- function(
    spacers, bsgenome,  method= c('Doench2014', 'Doench2016')[1],
    python = NULL, virtualenv = NULL, condaenv = NULL, 
    verbose = TRUE, plot = TRUE, 
    alpha_var = default_alpha_var(spacers)
){
    # Assert
    assert_is_all_of(spacers, 'GRanges')
    assert_is_a_string(method)
    assert_is_subset(method, c('Doench2014', 'Doench2016'))
    if (method %in% names(mcols(spacers))) mcols(spacers)[[method]] <- NULL

    # Add contextseq
    if (verbose)  cmessage('\tScore crispr spacers')
    spacers %<>% add_context(bsgenome, verbose = verbose)
    spacerdt  <- gr2dt(spacers)
    scoredt <- data.table(crisprcontext = unique(spacerdt$crisprcontext))
    
    # Score
    scorefun <- switch(method, Doench2014 = doench2014, Doench2016 = doench2016)
    scoredt[ , (method) := scorefun(scoredt$crisprcontext, verbose=verbose) ]

    # Merge back in
    mergedt  <- merge(spacerdt, scoredt, by='crisprcontext', sort=FALSE, all.x=TRUE)
    spacers <- dt2gr(mergedt, seqinfo = seqinfo(spacers))
    
    # Plot
    if (plot){
        scores   <- mcols(spacers)[[method]]
        tertiles <- quantile(scores, c(0.33, 0.66, 1))
        labels   <- sprintf('%s < %s (%s)', 
                        method, as.character(round(tertiles, 2)), names(tertiles))
        spacers$efficiency <- cut(scores, c(0, tertiles), labels)
        p <- plot_intervals(
                spacers, size_var = 'efficiency', alpha_var = alpha_var) + 
            ggplot2::scale_size_manual(values = c(0.1, 1, 2))
        print(p)
        spacers$efficiency <- NULL
    }
    
    # Return
    spacers
}

#' @rdname add_efficiency
#' @export
filter_efficient <- function(
    spacers, 
    bsgenome,  
    method= c('Doench2014', 'Doench2016')[1],
    cutoff,
    python     = NULL, 
    virtualenv = NULL, 
    condaenv   = NULL, 
    verbose    = TRUE, 
    plot       = TRUE,
    alpha_var  = default_alpha_var(gr)
){
    spacers %<>% add_efficiency(
                    bsgenome = bsgenome,  method = method, python = python, 
                    virtualenv = virtualenv, condaenv = condaenv, 
                    verbose = verbose, plot = plot, alpha_var = alpha_var)

    idx <- mcols(spacers)[[method]] > cutoff
    if (verbose){
        width <- nchar(length(idx))
        cmessage('\t\t%s ranges', formatC(length(idx), width = width))
        cmessage('\t\t%s ranges after filtering for %s > %s',
                 formatC(sum(idx), width = width), method, as.character(cutoff))
    }
    
    spacers %>% extract(idx)
}

default_alpha_var <- function(gr){
    if ('specific' %in% names(mcols(gr))) 'specific' else NULL
}