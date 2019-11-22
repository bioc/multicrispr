#============================================================================

#' Count target/genome matches
#' @param cas9seqs    character() or \code{\link[Biostrings]{XStringSet-class}}
#' @param targetseqs  character() or \code{\link[Biostrings]{XStringSet-class}}
#' @param bsgenome    \code{\link[BSgenome]{BSgenome-class}}
#' @param mismatch    number: number of allowed mismatches 
#' @param chromosomes character vector
#' @param verbose     logical(1)
#' @return numeric(length(cas9seqs))
#' @examples
#' # Read ranges and extend
#'     require(magrittr)
#'     bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'     targets <- bed_to_granges(bedfile, 'mm10')
#'     targets %<>% extend()
#' 
#' # Add seqs and find cas9s
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'     targets %<>% add_seq(bsgenome)
#'     cas9s <- find_crispr_sites(targets)
#' 
#' # Count matches
#'     cas9seqs <- cas9s$seq[1:10]
#'     count_target_matches(cas9seqs, targets$seq, 0)
#'     count_genome_matches(cas9seqs, bsgenome,    0, chromosomes = 'chrY')
#' @export
count_target_matches <- function(
    cas9seqs, 
    targetseqs, 
    mismatch, 
    verbose = TRUE
){
    
    # Assert
    assertive.types::assert_is_any_of(cas9seqs,   c('character', 'XStringSet'))
    assertive.types::assert_is_any_of(targetseqs, c('character', 'XStringSet'))
    assertive.types::assert_is_a_number(mismatch)
    assertive.sets::assert_is_subset(mismatch, c(0,1,2))
    assertive.types::assert_is_a_bool(verbose)
    
    # Count
    uniquecas9s <- unique(cas9seqs)
    starttime <- Sys.time()
    matches  <- rowSums(Biostrings::vcountPDict(
                            Biostrings::DNAStringSet(uniquecas9s),
                            Biostrings::DNAStringSet(targetseqs),
                            min.mismatch = mismatch,
                            max.mismatch = mismatch)) %>% 
                magrittr::set_names(uniquecas9s)
    
    # Return
    if (verbose) cmessage( '\t\t\tCount %d-mismatch hits in targets: %s',
                            mismatch,
                            format(signif(Sys.time() - starttime, 2)))
    unname(matches[cas9seqs])
}


#' @rdname count_target_matches
#' @export
count_genome_matches <- function(
    cas9seqs, 
    bsgenome, 
    mismatch,
    chromosomes = BSgenome::seqnames(bsgenome),
    verbose     = TRUE
){
    # Assert
    assertive.types::assert_is_any_of(cas9seqs, c('character', 'XStringSet'))
    assertive.types::assert_is_any_of(bsgenome, 'BSgenome')
    assertive.sets::assert_is_subset(mismatch, c(0,1,2))
    assertive.types::assert_is_character(chromosomes)
    assertive.types::assert_is_a_bool(verbose)

    # Comply
    . <- count <- NULL

    # Count
    starttime <- Sys.time()
    exclude  <- setdiff(BSgenome::seqnames(bsgenome), chromosomes)
    uniquecas9s <- unique(cas9seqs)
    matches  <- Biostrings::vcountPDict(
                    Biostrings::DNAStringSet(uniquecas9s),
                    bsgenome,
                    min.mismatch = mismatch,
                    max.mismatch = mismatch, 
                    exclude      = exclude, 
                    verbose      = verbose) %>% 
                data.table::as.data.table() %>% 
                magrittr::extract(, .(n = sum(count)), by ='index') %>%
                magrittr::extract2('n') %>% 
                magrittr::set_names(uniquecas9s)
    
    # Return
    if (verbose) cmessage( '\t\t\tCount %d-mismatch hits in genome      : %s',
                            mismatch, format(signif(Sys.time() - starttime, 2)))
    unname(matches[cas9seqs])
}


add_seqinfo <- function(gr, bsgenome){
    GenomeInfoDb::seqinfo(gr) <- GenomeInfoDb::seqinfo(bsgenome)
    gr
}


#' Filter for no offtargets
#' @param cas9s         \code{\link[GenomicRanges]{GRanges-class}}
#' @param targets       \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome      \code{\link[BSgenome]{BSgenome-class}}
#' @param mismatch      number: max number of mismatches to consider
#' @param offtargetchr  character vector: chromosomes for offtarget analysis
#' @param plot          logical(1)
#' @param verbose       logical(1)
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
#'     targets %<>% extract(GenomicRanges::seqnames(.)=='chrY')
#' 
#' # Add seqs and find cas9s
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'     targets %<>% add_seq(bsgenome)
#'     cas9s <- find_crispr_sites(targets)
#'    
#' # Filter for offtarget-free cas9s
#'     filter_no_offtargets(cas9s, targets, bsgenome, 0, offtargetchr = 'chrY')
#' @export
filter_no_offtargets <- function(cas9s, targets, bsgenome, mismatch = 2, 
    offtargetchr = GenomeInfoDb::standardChromosomes(targets), 
    plot = TRUE, verbose = TRUE){
    
    # Assert. Prepare
    assertive.types::assert_is_all_of(cas9s,    'GRanges')
    assertive.types::assert_is_all_of(targets,  'GRanges')
    assertive.types::assert_is_all_of(bsgenome, 'BSgenome')
    assertive.types::assert_is_a_number(mismatch)
    assertive.types::assert_is_character(offtargetchr)
    assertive.sets::assert_is_subset(offtargetchr, BSgenome::seqnames(bsgenome))
    assertive.sets::assert_is_subset('seq',names(GenomicRanges::mcols(cas9s)))
    assertive.sets::assert_is_subset('seq',names(GenomicRanges::mcols(targets)))
    targetseqs <- targets$seq
    cas9seqdt  <- data.table::data.table(seq = unique(cas9s$seq))
    
    # Count-store-filter for 0-2 mismatches
    if (verbose) message('\tFind cas9seq hits')
    for (mis in 0:mismatch){
        if (verbose) cmessage('\t\twith %d mismatch(es)', mis)
        # Count
        target_matches <- count_target_matches(
                            cas9seqdt$seq, targetseqs, mis, verbose = verbose)
        genome_matches <- count_genome_matches(cas9seqdt$seq, bsgenome, mis, 
                            chromosomes = offtargetchr, verbose = verbose)
        # Store
        cas9seqdt [ , (sprintf('matches%d', mis)) := target_matches ]
        
        # Filter
        assertive.base::assert_all_are_false(genome_matches < target_matches)
        idx <- genome_matches == target_matches
        if (verbose)   cmessage('\t\t\tKeep %d/%d offtarget-free cas9seqs', 
                                sum(idx), length(idx), mis)
        cas9seqdt %<>% magrittr::extract(idx)
    }
    
    # Filter for offtargetfree
    offtargetfree  <-   cas9seqdt %>% 
                        merge(data.table::as.data.table(cas9s), by = 'seq') %>% 
                        methods::as('GRanges') %>%  add_seqinfo(bsgenome)
    # Plot/Report
    if (plot)    plot_karyogram(GenomicRanges::GRangesList(
                                    target = targets, cas9 = cas9s, 
                                    offtargetfree = offtargetfree))
    if (verbose) cmessage('\tReturn %d cas9seqs across %d ranges',
                    length(unique(offtargetfree$seq)), length(offtargetfree))
    # Return
    offtargetfree
}
