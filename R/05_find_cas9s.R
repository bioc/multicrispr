
#' Convert GRanges into Sequences
#' @param granges \code{\link[GenomicRanges]{GRanges-class}}
#' @return character vector
#' @examples 
#' bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' granges <- read_bed(bedfile, 'mm10')
#' seqs(granges)
#' @export
seqs <- function(granges){
    assert_is_identical_to_true(is(granges, 'GRanges'))
    BSgenome::getSeq(get_bsgenome(granges), granges) %>%
    as.character()
}

#' Complement
#'
#' Adds inverse strand for each range
#'
#' @param gr \code{\link[GenomicRanges]{GRanges-class}}
#' @param plot     logical(1)
#' @param verbose  logical(1)
#' @return \code{\link[GenomicRanges]{GRanges-class}}, twice as long as input
#' @examples
#' require(magrittr)
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' gr <- read_bed(bedfile, 'mm10', plot = FALSE)
#' complement(gr)
#' @export
complement <- function(gr, plot = TRUE, verbose = TRUE){
    complements <- invertStrand(gr)
    newranges <- c(gr, complements)
    txt <- sprintf('\t\t%d ranges after adding inverse strands',
                   length(newranges))
    if (plot){
        plot_intervals(
            GRangesList(original = gr, complements = complements),
            title = txt)
    }
    if (verbose) cmessage(txt)
    newranges
}




#=============================================================================
# Find cas9 ranges
#=============================================================================

#' Find cas9 sites in targetranges
#' @param gr        \code{\link[GenomicRanges]{GRanges-class}}
#' @param inclcompl logical(1): include complementary strands in search?
#' @param verbose   logical(1): report verbosely?
#' @return   \code{\link[GenomicRanges]{GRanges-class}}
#' @examples
#' bedfile  <- system.file('extdata/SRF.bed', package='multicrispr')
#' targetranges <- slop_fourways(read_bed(bedfile, 'mm10'))
#' find_cas9s(targetranges)
#' @export 
find_cas9s <- function(gr, inclcompl = TRUE, verbose = TRUE){

    # Assert
    assert_is_identical_to_true(is(gr, 'GRanges'))
    assert_is_a_bool(verbose)
    
    # Add complementary strands
    if (verbose) message('\tFind N{20}NGG cas9seqs')
    if (inclcompl) gr %<>% complement(plot = FALSE, verbose = verbose)
    
    # Comply
    start <- substart <- cas9start <- NULL
    end <- subend <- cas9end <- strand <- seqnames <- NULL
    
    # Find cas9s in targetranges
    targetdt <- data.table::as.data.table(gr)
    targetdt [ , seqs := seqs(gr) ]
    res <- targetdt$seqs %>% stringi::stri_locate_all_regex('[ACGT]{21}GG')
    cextract1 <- function(y) y[, 1] %>% paste0(collapse=';')
    cextract2 <- function(y) y[, 2] %>% paste0(collapse=';')
    targetdt [ , substart := vapply( res, cextract1, character(1)) ]
    targetdt [ , subend   := vapply( res, cextract2, character(1)) ]
    
    # Rm cas9-free targetranges
    idx <- targetdt[, substart == 'NA']
    if (sum(idx)>0){
        if (verbose)  cmessage('\t\tRm %d ranges with no cas9sites', sum(idx)) 
        targetdt %<>% extract(!idx)
    }

    # Transform into cas9ranges
    cas9ranges <- tidyr::separate_rows(targetdt, substart, subend) %>%
        data.table::data.table()                                   %>% 
        extract(, substart := as.numeric(substart))                %>% 
        extract(, subend   := as.numeric(subend))                  %>% 
        extract(, seqs     := substr(seqs, substart, subend))      %>%
        extract( strand=='+', cas9start := start + substart - 1  ) %>% 
        extract( strand=='+', cas9end   := start + subend   - 1  ) %>%
        extract( strand=='-', cas9start := end   - subend   + 1  ) %>%
        extract( strand=='-', cas9end   := end   - substart + 1  ) %>%
        extract(, list(seqnames = seqnames, start  = cas9start, end = cas9end,  
                       strand   = strand,  seqs    = seqs) ) %>% 
        unique() %>% as('GRanges') %>% add_seqinfo(get_bsgenome(gr))
    
    # Return
    if (verbose)   cmessage('\t\t%d cas9 seqs across %d ranges', 
                            length(unique(seqs(cas9ranges))), 
                            length(cas9ranges))
    return(cas9ranges)
}

#============================================================================

#' Count target matches
#' @param cas9seqs   character() or \code{\link[Biostrings]{XStringSet-class}}
#' @param targetseqs character() or \code{\link[Biostrings]{XStringSet-class}}
#' @param mismatch   number: number of allowed mismatches 
#' @param verbose    logical(1)
#' @return numeric(length(cas9seqs))
#' @examples
#' 
#' # Read target ranges
#' require(magrittr)
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' targetranges <- read_bed(bedfile, 'mm10')
#' targetranges %<>% slop_fourways %>% extract(1:10)
#' 
#' # Find cas9ranges
#' cas9ranges <- find_cas9s(targetranges)
#' 
#' # Count target matches
#' bsgenome <- get_bsgenome(cas9ranges)
#' count_target_matches(
#'     cas9seqs   = getSeq(bsgenome, cas9ranges),
#'     targetseqs = getSeq(bsgenome, targetranges),
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
    assert_is_any_of(cas9seqs,   c('character', 'XStringSet'))
    assert_is_any_of(targetseqs, c('character', 'XStringSet'))
    assert_is_a_number(mismatch)
    assert_is_subset(mismatch, c(0,1,2))
    assert_is_a_bool(verbose)
    
    # Count
    starttime <- Sys.time()
    matches  <- rowSums(vcountPDict(DNAStringSet(cas9seqs),
                                    DNAStringSet(targetseqs),
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
#' targetranges <- read_bed(bedfile, 'mm10')
#' targetranges %<>% slop_fourways %>% extract(1:10)
#' 
#' # Find cas9ranges
#' cas9ranges <- find_cas9s(targetranges)
#' 
#' # Count genome matches
#' bsgenome <- get_bsgenome(cas9ranges)
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
    chromosomes = seqnames(bsgenome),
    verbose     = TRUE
){
    # Assert
    assert_is_any_of(cas9seqs, c('character', 'XStringSet'))
    assert_is_any_of(bsgenome, 'BSgenome')
    assert_is_subset(mismatch, c(0,1,2))
    assert_is_character(chromosomes)
    assert_is_a_bool(verbose)

    # Comply
    . <- count <- NULL

    # Count
    starttime <- Sys.time()
    exclude  <- setdiff(seqnames(bsgenome), chromosomes)
    matches  <- vcountPDict(DNAStringSet(cas9seqs),
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

#' Get (canonical) chromosomes
#' @param gr \code{\link[GenomicRanges]{GRanges-class}}
#' @return character vector
#' @examples 
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' gr <- read_bed(bedfile, 'mm10')
#' chromosomes(gr)
#' canonicalchr(gr)
#' @export
chromosomes <- function(gr){
    seqnames(seqinfo(gr))
}

#' @rdname chromosomes
#' @export
canonicalchr <- function(gr){
    chromosomes(gr) %>%
    extract(stringi::stri_detect_fixed(., '_random', negate = TRUE)) %>% 
    extract(stringi::stri_detect_fixed(., 'chrUn', negate = TRUE))
}

#' Find offtarget-free cas9 sites in targetranges
#' @param targetranges  \code{\link[GenomicRanges]{GRanges-class}}
#' @param mismatch        number: max number of mismatches to consider
#' @param offtargetchr character vector: chromosomes for offtarget analysis, 
#'                     probably generated with chromosomes(targetranges) or 
#'                     canonicalchr(targetranges)
#' @param verbose  logical(1)
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#'         mcols(GRanges) contains sequences and match counts:
#'             matches0 = perfect match counts
#'             matches1 = single mismatch counts
#'             matches2 = double mismatch counts
#' @examples
#' # Note: restricting example to 'chrY' only to keep it fast
#' bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' targetranges <- read_bed(bedfile, 'mm10', plot = FALSE) %>% 
#'                 extract(seqnames(.)=='chrY') %>% 
#'                 slop_fourways()
#' find_offtargetfree_cas9s(targetranges, 0, offtargetchr = 'chrY')
#' @export
find_offtargetfree_cas9s <- function(
    targetranges, 
    mismatch     = 2,
    offtargetchr = canonicalchr(targetranges),
    verbose      = TRUE
){
    # Assert
    assertive.types::assert_is_a_number(mismatch)

    # Prepare
    bsgenome   <- get_bsgenome(targetranges)
    targetseqs <- seqs(targetranges)
    cas9ranges <- targetranges %>% find_cas9s(verbose = verbose)
    cas9seqdt  <- data.table::data.table(seqs = unique(cas9ranges$seqs))
    
    # Count-store-filter for 0-2 mismatches
    if (verbose) message('\tFind cas9seq hits')
    for (mis in 0:mismatch){
        if (verbose) cmessage('\t\twith %d mismatch(es)', mis)
        # Count
        target_matches  <-  count_target_matches(
                                cas9seqdt$seqs, targetseqs, mis, 
                                verbose = verbose)
        genome_matches  <-  count_genome_matches(
                                cas9seqdt$seqs, bsgenome, mis, 
                                chromosomes = offtargetchr, 
                                verbose = verbose)
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
    specificranges  <-  cas9seqdt %>% 
                        merge(data.table::as.data.table(cas9ranges), 
                              by = 'seqs') %>% 
                        as('GRanges') %>% 
                        add_seqinfo(bsgenome)
        
    if (verbose)    cmessage('Returning %d cas9seqs across %d ranges',
                            length(unique(specificranges$seqs)), 
                            length(specificranges))
    specificranges
        
}
