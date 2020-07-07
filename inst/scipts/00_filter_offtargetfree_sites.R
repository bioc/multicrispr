#============================================================================

#' Count target/genome matches
#' @param spacerseqs  character vector or 
#'                    \code{\link[Biostrings]{XStringSet-class}}
#' @param targetseqs  character vector or 
#'                    \code{\link[Biostrings]{XStringSet-class}}
#' @param pam         string (default 'NGG'): sequence following spacerseqs
#' @param bsgenome    \code{\link[BSgenome]{BSgenome-class}}
#' @param mismatch    number: number of allowed mismatches 
#' @param include     character vector: offtarget analysis chromosomes
#' @param verbose     TRUE (default) or FALSE
#' @return numeric(length(spacerseqs))
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
#'     spacerseqs <- spacers$spacer
#'     # count_genome_matches(spacerseqs, bsgenome, mismatch=0) # 6 mins
#'     # count_genome_matches(spacerseqs, bsgenome, mismatch=1) # 32 mins
#'     # count_genome_matches(spacerseqs, bsgenome, mismatch=2)
#'     # count_genome_matches(spacerseqs, bsgenome, mismatch=3)
#'     
#' # TFBS example
#' #-------------
#'     require(magrittr)
#'     bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'     targets <- extend(bed_to_granges(bedfile, 'mm10'))
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'     spacers <- find_spacers(targets, bsgenome)
#'     spacerseqs <- spacers$spacer
#'     count_genome_matches(spacerseqs, bsgenome, mismatch = 0, include = 'chrY')
#' 
#' # Count matches
#'     spacerseqs <- sites$seq[1:10]
#'     targetseqs <- targets$seq
#'     count_target_matches(spacerseqs, targetseqs, 0)
#'     count_genome_matches(spacerseqs, bsgenome,    0, include = 'chrY')
#' @export
count_target_matches <- function(
    spacerseqs, 
    targetseqs, 
    mismatch, 
    verbose = TRUE
){
    
    # Assert
    assert_is_any_of(spacerseqs, c('character', 'XStringSet'))
    assert_is_any_of(targetseqs, c('character', 'XStringSet'))
    assert_is_a_number(mismatch)
    assert_is_subset(mismatch, 0:4)
    assert_is_a_bool(verbose)
    
    # Count
    unique_crisprseqs <- unique(spacerseqs)
    starttime <- Sys.time()
    matches  <- rowSums(Biostrings::vcountPDict(
                            Biostrings::DNAStringSet(unique_crisprseqs),
                            Biostrings::DNAStringSet(targetseqs),
                            min.mismatch = mismatch,
                            max.mismatch = mismatch)) %>% 
                set_names(unique_crisprseqs)
    
    # Return
    if (verbose) cmessage( '\t\t\tCount %d-mismatch hits in targets: %s',
                            mismatch,
                            format(signif(Sys.time() - starttime, 2)))
    unname(matches[spacerseqs])
}


#' Convert "included" into "excluded" seqnames
#'
#' @details Convert "included" into "excluded" seqnames. 
#'          Do this in regex format as required by 
#' \code{\link[BSgenome]{BSgenome-utils}}) \code{vcountPDict} and 
#' \code{vmatchPDict}
#' @param bsgenome \code{\link[BSgenome]{BSgenome-class}}
#' @param include character vector with names of included chromosomes
#' @return character vector
#' @examples 
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#' include_to_exclude(bsgenome)
#' @export
include_to_exclude <- function(
    bsgenome, 
    include = GenomeInfoDb::standardChromosomes(bsgenome)
){
    assert_is_all_of(bsgenome, 'BSgenome')
    assert_is_character(include)
    
    # needs to be a regex: https://support.bioconductor.org/p/126593/
    sprintf('^%s$', setdiff(seqnames(bsgenome), include))
}

    explode <- function(x) unlist(strsplit(x, character(0)))
    paste_dtcols <- function(x) x [, do.call(paste0, .SD) ]
    expand_iupac_ambiguities <- function(x){
        paste_dtcols( 
            as.data.table( 
                lapply(Biostrings::IUPAC_CODE_MAP[explode(x)], explode)))
    }

    vsprintf <- function(fmt, ..., first_slowest = TRUE){
       . <- NULL
       arguments <- list(...)
       n <- arguments %>% vapply(length, integer(1)) %>% Reduce(max, .)
       argdf <- arguments %>%
               (function(x) if (first_slowest) rev(x) else x) %>%
                expand.grid() %>%
               (function(x) if (first_slowest) x %>% magrittr::extract(, ncol(.):1) else x)
       fun <- function(...) sprintf(fmt, ...)
       do.call(fun, argdf)
    }
    

.match_genome <- function(seqs, bsgenome, exclude, mismatch){
        res <- vcountPDict(seqs, bsgenome, min.mismatch = mismatch, 
                                max.mismatch = mismatch, exclude = exclude)
        res %<>% data.table::as.data.table()
        counts = res[, .(N = sum(count)), by = index]$N
}


.match_genome_slow <- function(seqs, bsgenome, exclude, mismatch){
    # Parallel R-loop equivalent to .match_genome (serial C loops)
    # Is slower, more complicated, and more resource-consuming.
    # So no longer used.
    unlist(BiocParallel::bplapply(
        seqs, 
        function(x, bsgenome, exclude, mismatch){
            sum(BSgenome::vcountPattern(
                    pattern      = Biostrings::DNAString(x), 
                    subject      = bsgenome, 
                    min.mismatch = mismatch, 
                    max.mismatch = mismatch, 
                    exclude      = exclude)$count)},
        bsgenome = bsgenome, exclude = exclude, mismatch = mismatch))
}


#' @rdname count_target_matches
#' @export
count_genome_matches <- function(
    spacerseqs, 
    bsgenome, 
    mismatch = 0,
    pam      = 'NGG',
    include  = GenomeInfoDb::standardChromosomes(bsgenome),
    verbose  = TRUE, 
    joint    = FALSE
){
    # Assert
    assert_is_any_of(spacerseqs, c('character', 'XStringSet'))
    assert_is_a_string(pam)
    assert_is_any_of(bsgenome, 'BSgenome')
    assert_is_subset(mismatch, 0:4)
    assert_is_character(include)
    assert_is_a_bool(verbose)

    # Comply
    . <- count <- NULL

    # Count
    pamseqs <- expand_iupac_ambiguities(pam)
    sitedt <- data.table(spacer = rep(spacerseqs, each = length(pamseqs)), 
                        pam     = rep(pamseqs, times = length(spacerseqs))) %>% 
              extract(, crispr := paste0(spacer, pamseqs) )
    seqdt <- data.table(crispr = unique(sitedt$crispr))
    
    starttime <- Sys.time()
    excl <- include_to_exclude(bsgenome, include)
    seqdt[ , nmatches := .match_genome(crispr, bsgenome, excl, mismatch) ]
    
    # Return
    if (verbose) cmessage( '\t\t\tCount %d-mismatch hits in genome      : %s',
                            mismatch, format(signif(Sys.time() - starttime, 2)))
    sitedt %>%  merge(seqdt, by = 'crispr') %>% 
                extract(, .(n = sum(nmatches)), by = 'spacer') %>% 
                extract(spacerseqs, n, on = 'spacer')
}


add_seqinfo <- function(gr, bsgenome){
    seqinfo(gr) <- seqinfo(bsgenome)
    gr
}


#' Filter for offtarget-free crispr sites
#' @param sites    \code{\link[GenomicRanges]{GRanges-class}}: crispr sites
#' @param targets  \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome \code{\link[BSgenome]{BSgenome-class}}
#' @param mismatch number: max number of mismatches to consider
#' @param include  character vector: offtarget analysis chromosomes
#' @param plot     TRUE (default) or FALSE
#' @param verbose  TRUE (default) or FALSE
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#'         mcols(GRanges) contains sequences and match counts:
#'             matches0 = perfect match counts
#'             matches1 = single mismatch counts
#'             matches2 = double mismatch counts
#' @examples
#' # Read ranges, extend, and restrict to 'chrY' (to keep example fast)
#'     require(magrittr)
#'     bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'     targets <- bed_to_granges(bedfile, 'mm10')
#'     targets %<>% extend()
#'     targets %<>% extract(seqnames(.)=='chrY')
#' 
#' # Add seqs and find crispr sites
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'     targets %<>% add_seq(bsgenome)
#'     sites <- find_spacers(targets, bsgenome)
#'    
#' # Filter for offtarget-free sites
#'     filter_offtargetfree_sites(
#'         sites, targets, bsgenome, 0, include = 'chrY')
#' @export
filter_offtargetfree_sites <- function(sites, targets, bsgenome, mismatch = 2, 
    include = standardChromosomes(targets), plot = TRUE, verbose = TRUE){
    
    # Assert. Prepare
    assert_is_all_of(sites,    'GRanges')
    assert_is_all_of(targets,  'GRanges')
    assert_is_all_of(bsgenome, 'BSgenome')
    assert_is_a_number(mismatch)
    assert_is_character(include)
    assert_is_subset(include, seqnames(bsgenome))
    assert_is_subset('seq',names(mcols(sites)))
    assert_is_subset('seq',names(mcols(targets)))
    targetseqs <- targets$seq
    crisprseq_dt  <- data.table(seq = unique(sites$seq))
    
    # Count-store-filter for 0-2 mismatches
    if (verbose) message('\tFind crispr seq (mis)matches')
    for (mis in 0:mismatch){
        if (verbose) cmessage('\t\twith %d mismatch(es)', mis)
        # Count
        target_matches <- count_target_matches(crisprseq_dt$seq, 
                                            targetseqs, mis, verbose = verbose)
        genome_matches <- count_genome_matches(crisprseq_dt$seq, bsgenome, mis, 
                            include = include, verbose = verbose)
        # Store
        crisprseq_dt [ , (sprintf('matches%d', mis)) := target_matches ]
        
        # Filter
        assert_all_are_false(genome_matches < target_matches)
        idx <- genome_matches == target_matches
        if (verbose)   cmessage('\t\t\tKeep %d/%d offtarget-free crispr seqs', 
                                sum(idx), length(idx), mis)
        crisprseq_dt %<>% extract(idx)
    }
    
    # Filter for offtargetfree sites
    offtargetfree  <-   crisprseq_dt %>% 
                        merge(as.data.table(sites), by = 'seq') %>% 
                        as('GRanges') %>%  add_seqinfo(bsgenome)
    # Plot/Report
    if (plot)    plot_karyogram(GenomicRanges::GRangesList(
                                    target = targets, cas9 = sites, 
                                    offtargetfree = offtargetfree))
    if (verbose) cmessage('\tReturn %d crispr seqs across %d ranges',
                    length(unique(offtargetfree$seq)), length(offtargetfree))
    # Return
    offtargetfree
}
