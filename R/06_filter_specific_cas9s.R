#============================================================================

#' Count target matches
#' @param cas9seqs   character() or \code{\link[Biostrings]{XStringSet-class}}
#' @param targetseqs character() or \code{\link[Biostrings]{XStringSet-class}}
#' @param mismatch   number: number of allowed mismatches 
#' @param verbose    logical(1)
#' @return numeric(length(cas9seqs))
#' @examples
#' # Read target ranges
#' require(magrittr)
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' targetranges <- bed_to_granges(bedfile, bsgenome)
#' targetranges %<>% extend() %>% extract(1:10)
#' 
#' # Find cas9ranges
#' cas9ranges <- find_cas9s(targetranges)
#' 
#' # Count target matches
#' count_target_matches(
#'     cas9seqs   = seqs(cas9ranges, bsgenome),
#'     targetseqs = seqs(targetranges, bsgenome),
#'     mismatch   = 0, 
#'     verbose    = TRUE)
#' @seealso \code{\link{count_target_matches}}, \code{\link{vcountPDict}}
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
    starttime <- Sys.time()
    matches  <- rowSums(Biostrings::vcountPDict(
                            Biostrings::DNAStringSet(cas9seqs),
                            Biostrings::DNAStringSet(targetseqs),
                            min.mismatch = mismatch,
                            max.mismatch = mismatch))
    # Return
    if (verbose) cmessage( '\t\t\tCount %d-mismatch hits in targetranges: %s',
                            mismatch,
                            format(signif(Sys.time() - starttime, 2)))
    matches
}


#' Count genome matches
#' @param cas9seqs     character() or \code{\link[Biostrings]{XStringSet-class}}
#' @param bsgenome     \code{\link[BSgenome]{BSgenome-class}}
#' @param mismatch     number: number of allowed mismatches 
#' @param chromosomes  character vector
#' @param verbose      logical(1)
#' @return numeric(length(cas9seqs))
#' @examples
#' 
#' # Read target ranges
#' require(magrittr)
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' targetranges <- bed_to_granges(bedfile, bsgenome)
#' targetranges %<>% extend() %>% extract(1:10)
#' 
#' # Find cas9ranges
#' cas9ranges <- find_cas9s(targetranges)
#' 
#' # Count genome matches
#' count_genome_matches(
#'     cas9seqs    = getSeq(bsgenome, cas9ranges),
#'     bsgenome    = bsgenome, 
#'     mismatch    = 0, 
#'     chromosomes = 'chrY',
#'     verbose     = TRUE)
#' @seealso \code{\link{count_genome_matches}}, \code{\link{BSgenome-utils}}
#' @export
count_genome_matches <- function(
    cas9seqs, 
    bsgenome, 
    mismatch,
    chromosomes = GenomeInfoDb::seqnames(bsgenome),
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
    exclude  <- setdiff(GenomeInfoDb::seqnames(bsgenome), chromosomes)
    matches  <- Biostrings::vcountPDict(
                    Biostrings::DNAStringSet(cas9seqs),
                    bsgenome,
                    min.mismatch = mismatch,
                    max.mismatch = mismatch, 
                    exclude      = exclude) %>% 
                data.table::as.data.table() %>% 
                extract(, .(n = sum(count)), by ='index') %>%
                extract2('n') 

    # Return
    if (verbose) cmessage( '\t\t\tCount %d-mismatch hits in genome      : %s',
                            mismatch, format(signif(Sys.time() - starttime, 2)))
    matches
}

#' Filter for specific (i.e. offtarget-free) cas9 sites
#' @param cas9ranges    \code{\link[GenomicRanges]{GRanges-class}}
#' @param targetranges  \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome      \code{\link[BSgenome]{BSgenome-class}}
#' @param mismatch        number: max number of mismatches to consider
#' @param offtargetchr character vector: chromosomes for offtarget analysis, 
#'          probably generated with genomeInfoDb::seqlevels(targetranges) or 
#'          canonicalseqlevels(targetranges)
#' @param verbose  logical(1)
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#'         mcols(GRanges) contains sequences and match counts:
#'             matches0 = perfect match counts
#'             matches1 = single mismatch counts
#'             matches2 = double mismatch counts
#' @examples
#' # Note: restricting example to 'chrY' only to keep it fast
#' require(magrittr)
#' bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' targetranges <- bed_to_granges(bedfile, bsgenome, plot = FALSE) %>% 
#'                 extract(GenomeInfoDb::seqnames(.)=='chrY') %>% 
#'                 extend()
#' cas9ranges <- find_cas9s(targetranges, bsgenome)
#' find_offtargetfree_cas9s(cas9ranges, targetranges, bsgenome, mismatch = 0, 
#'                          offtargetchr = 'chrY')
#' @export
filter_no_offtargets <- function(
    cas9ranges, 
    targetranges,
    bsgenome, 
    mismatch     = 2,
    offtargetchr = canonicalseqlevels(targetranges),
    verbose      = TRUE
){
    # Assert
    assertive.types::assert_is_a_number(mismatch)

    # Prepare
    targetseqs <- seqs(targetranges)
    cas9seqdt  <- data.table::data.table(seqs = unique(cas9ranges$seqs))
    
    # Count-store-filter for 0-2 mismatches
    if (verbose) message('\tFind cas9seq hits')
    for (mis in 0:mismatch){
        if (verbose) cmessage('\t\twith %d mismatch(es)', mis)
        # Count
        target_matches <- count_target_matches(
                            cas9seqdt$seqs, targetseqs, mis, verbose = verbose)
        genome_matches <- count_genome_matches(
                            cas9seqdt$seqs, bsgenome, mis, 
                            chromosomes = offtargetchr, verbose = verbose)
        # Store
        cas9seqdt [ , (sprintf('matches%d', mis)) := target_matches ]
        
        # Filter
        assertive.base::assert_all_are_false(genome_matches < target_matches)
        idx <- genome_matches == target_matches
        if (verbose)   cmessage('\t\t\tKeep %d/%d target-specific cas9seqs', 
                                sum(idx), length(idx), mis)
        cas9seqdt %<>% extract(idx)
    }
    
    # Return
    specificranges  <-  merge(cas9seqdt, 
                            data.table::as.data.table(cas9ranges), 
                            by = 'seqs') %>% 
                        as('GRanges') %>% 
                        add_seqinfo(bsgenome)
        
    if (verbose)    cmessage('Returning %d cas9seqs across %d ranges',
                            length(unique(specificranges$seqs)), 
                            length(specificranges))
    specificranges
}
