
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
#' @noRd
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

# Note: mclapply vs. bplapply
#     I first used BiocParallel::bplapply().
#     That works great on linux, but it fails on windows.
#     support.bioconductor.org/p/92587/ explains likely why
#     Windows doensn't allow forking (where a child process inherits parent env)
#     Instead bplapply defaults to multithreading.
#     In that scenario, the reticulation env is not correctly passed.
#    mclapply seems like best choice: fork on linux - execute serially on win
doench2016 <- function(
    contextseqs, 
    chunksize = 10000,
    verbose   = TRUE
){
    # Assert
    is_identical_to_true(reticulate::py_module_available('azimuth'))
    assert_is_character(contextseqs)
    assert_all_are_true(nchar(contextseqs)==30)
    assert_is_a_bool(verbose)
    
    # Message
    if (verbose)  message('\t\tScore contextseqs with Doench2016 (azimuth)')
    start_time <- Sys.time()
    
    # Score
    azi <- reticulate::import('azimuth.model_comparison', delay_load = TRUE)
    nchunks <- ceiling(length(contextseqs) / chunksize)
    contextchunks <- split(
                        contextseqs, ceiling(seq_along(contextseqs)/chunksize))
    txt <- paste0('\t\tRun Doench2016 %d times on %d-seq chunks ', 
                'to preserve memory')
    cmessage(txt, length(contextchunks), chunksize)
    mc.cores <- if (is_windows()) 1 else max(1, parallel::detectCores()-2)
    doench2016scores <- unlist(parallel::mclapply(contextchunks, 
            function(x){
                reticulate::py_suppress_warnings(
                    azi$predict( reticulate::np_array(x), 
                            aa_cut                 = NULL, 
                            percent_peptide        = NULL, 
                            model                  = NULL, 
                            model_file             = NULL, 
                            pam_audit              = TRUE, 
                            length_audit           = TRUE, 
                            learn_options_override = NULL))}, 
            mc.cores = mc.cores))
    
    # Return
    end_time <- Sys.time()
    if (verbose) cmessage('\t\tCompleted in %s', 
                        format(end_time - start_time, digits = 2))
    doench2016scores
} 


#' Add on-target efficiency scores
#' 
#' Add Doench2014 or Doench2016 on-target efficiency scores
#' 
#' \code{add_ontargets} adds efficiency scores
#' \code{filter_ontargets} adds efficiency scores and filters on them
#' 
#' @param spacers  \code{\link[GenomicRanges]{GRanges-class}}: spacers
#' @param bsgenome \code{\link[BSgenome]{BSgenome-class}}
#' @param ontargetmethod   'Doench2014' (default) or 'Doench2016'
#'                 (requires non-NULL argument python, virtualenv, or condaenv)
#' @param chunksize Doench2016 is executed in chunks of chunksize
#' @param verbose   TRUE (default) or FALSE
#' @param plot      TRUE (default) or FALSE
#' @param ...       passed to \code{\link{plot_intervals}}
#' @return numeric vector
#' @examples
#' # Install azimuth 
#' #----------------
#'     ## With reticulate
#'     # require(reticulate)
#'     # conda_create('azienv', c('python=2.7'))
#'     # use_condaenv('azienv')
#'     # py_install(c('azimuth', 'scikit-learn==0.17.1', 'biopython=='1.76'), 
#'     #            'azienv', pip = TRUE)
#'     
#'     ## Directly
#'     # conda create --name azienv python=2.7
#'     # conda activate azienv
#'     # pip install scikit-learn==0.17.1
#'     # pip install biopython==1.76
#'     # pip install azimuth
#'     
#' # PE example
#' #-----------
#'     require(magrittr)
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#'     targets <- char_to_granges(c(PRNP = 'chr20:4699600:+',             # snp
#'                                  HBB  = 'chr11:5227002:-',             # snp
#'                                  HEXA = 'chr15:72346580-72346583:-',   # del
#'                                  CFTR = 'chr7:117559593-117559595:+'), # ins
#'                                bsgenome)
#'     spacers <- find_primespacers(targets, bsgenome, ontargetmethod=NULL, 
#'                                 offtargetmethod=NULL)
#'     spacers %<>% score_ontargets(bsgenome, 'Doench2014')
#'     # reticulate::use_condaenv('azienv')
#'     # reticulate::import('azimuth')
#'     # spacers %<>% score_ontargets(bsgenome, 'Doench2016')
#'     
#' # TFBS example
#' #-------------
#'     bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'     targets <- extend(bed_to_granges(bedfile, 'mm10'))
#'     spacers <- find_spacers(targets, bsgenome, ontargetmethod=NULL, 
#'                             offtargetmethod=NULL)
#'     spacers %<>% score_ontargets(bsgenome, 'Doench2014')
#'     # reticulate::use_condaenv('azienv')
#'     # reticulate::import('azimuth')
#'     # spacers %>%  score_ontargets(bsgenome, 'Doench2016')
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
score_ontargets <- function(
    spacers, bsgenome,  ontargetmethod= c('Doench2014', 'Doench2016')[1],
    chunksize = 10000, verbose = TRUE, plot = TRUE, ...
){
    # Assert
    crisprcontext <- . <- NULL
    assert_is_all_of(spacers, 'GRanges')
    if (is.null(ontargetmethod))  return(spacers)
    assert_is_a_string(ontargetmethod)
    assert_is_subset(ontargetmethod, c('Doench2014', 'Doench2016'))
    if (ontargetmethod %in% names(mcols(spacers))){
        mcols(spacers)[[ontargetmethod]] <- NULL}

    # Add contextseq
    if (verbose)  cmessage('\tScore ontargets')
    spacers %<>% add_context(bsgenome, verbose = verbose)
    spacers %<>% extract(!stri_detect_fixed(.$crisprcontext, 'N'))
    spacerdt  <- gr2dt(spacers)  # context includes spacer + pam already
    scoredt <- data.table(crisprcontext = unique(spacerdt$crisprcontext))
    
    # Score
    scores <- switch(ontargetmethod, 
        Doench2014 = doench2014(scoredt$crisprcontext, verbose=verbose), 
        Doench2016 = doench2016(scoredt$crisprcontext, chunksize=chunksize, 
                                verbose=verbose))
    scoredt[ , (ontargetmethod) := scores ]

    # Merge back in
    mergedt  <- merge(spacerdt, scoredt,
                    by='crisprcontext', sort=FALSE, all.x=TRUE)
    mergedt[, crisprcontext := NULL]
    spacers <- dt2gr(mergedt, seqinfo = seqinfo(spacers))
    
    # Plot
    if (plot)   print(plot_intervals(spacers, ...))
    
    # Return
    spacers
}


default_alpha_var <- function(gr){
    if ('off' %in% names(mcols(gr))) 'off' else NULL
}
