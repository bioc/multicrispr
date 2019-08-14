#==============================================================================
# Get sequences
#==============================================================================

GRanges      <- methods::getClassDef('GRanges',      package = 'GenomicRanges')
DNAStringSet <- methods::getClassDef('DNAStringSet', package = 'Biostrings')


#' Get/set sequence values
#' @param granges GenomicRanges::GRanges
#' @param value character vector or DNAStringSet
#' @return DNAStringSet (get) or GRanges (set)
#' @examples 
#' bedfile  <- system.file('extdata/SRF_sites.bed', package = 'multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' granges <- read_bed(bedfile, bsgenome)
#' seqs(granges)
#' @export
seqs <- function(granges){
    assertive.base::assert_is_identical_to_true(is(granges, 'GRanges'))
    BSgenome::getSeq(get_bsgenome(granges), granges) %>% 
    as.character()
}


#=============================================================================
# Find cas9 ranges
#=============================================================================

#' Find cas9 ranges in targetranges
#' @param targetranges  GenomicRanges::GRanges
#' @param specific      logical(1)
#' @param verbose       logical(1)
#' @return GenomicRanges::GRanges
#' @examples
#' require(magrittr)
#' bedfile  <- system.file('extdata/SRF_sites.bed', package='multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' targetranges <- read_bed(bedfile, bsgenome) %>% flank_fourways()
#' find_cas9ranges(targetranges)
#' @export 
find_cas9ranges <- function(
    targetranges, 
    verbose = TRUE
){
    # Assert
    assertive.base::assert_is_identical_to_true(is(targetranges, 'GRanges'))
    assertive.types::assert_is_a_bool(verbose)
    
    # Find cas9s in targetranges
    if (verbose) message('\tFind N{20}NGG cas9seqs')
    targetdt <- data.table::as.data.table(targetranges)
    targetdt [ , seqs := seqs(targetranges) ]
    res <- targetdt$seqs %>% stringi::stri_locate_all_regex('[ACGT]{21}GG')
    cextract1 <- function(y) y[, 1] %>% paste0(collapse=';')
    cextract2 <- function(y) y[, 2] %>% paste0(collapse=';')
    targetdt [ , substart := vapply( res, cextract1, character(1)) ]
    targetdt [ , subend   := vapply( res, cextract2, character(1)) ]
    
    # Rm cas9-free targetranges
    idx <- targetdt[, substart == 'NA']
    if (sum(idx)>0){
        if (verbose)  cmessage('\t\tRm %d targetranges with no cas9sites', sum(idx)) 
        targetdt %<>% extract(!idx)
    }

    # Transform into cas9ranges
    cas9ranges <- targetdt                                                 %>% 
                tidyr::separate_rows(substart, subend)                     %>%
                data.table::data.table()                                   %>% 
                extract(, substart := as.numeric(substart))                %>% 
                extract(, subend   := as.numeric(subend))                  %>% 
                extract(, seqs     := substr(seqs, substart, subend))      %>%
                extract( strand=='+', cas9start := start + substart - 1  ) %>% 
                extract( strand=='+', cas9end   := start + subend   - 1  ) %>%
                extract( strand=='-', cas9start := end   - subend   + 1  ) %>%
                extract( strand=='-', cas9end   := end   - substart + 1  ) %>%
                extract( , list(seqnames = seqnames, 
                                start    = cas9start, 
                                end      = cas9end, 
                                strand   = strand, 
                                seqs     = seqs) ) %>% 
                unique() %>% 
                as.granges(get_bsgenome(targetranges))
    # Return
    if (verbose)   cmessage('\t\t%d cas9 seqs across %d targetranges', 
                            length(unique(seqs(cas9ranges))), 
                            length(cas9ranges))
    return(cas9ranges)
}

#============================================================================

count_target_matches <- function(cas9seqs, targetseqs, mismatch, verbose){
    
    starttime <- Sys.time()
    matches  <- Biostrings::vcountPDict(
                    Biostrings::DNAStringSet(cas9seqs),
                    Biostrings::DNAStringSet(targetseqs),
                    min.mismatch = mismatch,
                    max.mismatch = mismatch) %>% 
                rowSums()
    if (verbose) cmessage( '\t\t\tCount %d-mismatch hits in targetranges: %s',
                            mismatch,
                            format(signif(Sys.time() - starttime, 2)))
    matches
}


count_genome_matches <- function(cas9seqs, bsgenome, mismatch, verbose){
    
    starttime <- Sys.time()
    matches   <-    Biostrings::vcountPDict(
                        Biostrings::DNAStringSet(cas9seqs),
                        bsgenome,
                        min.mismatch = mismatch,
                        max.mismatch = mismatch) %>% 
                    data.table::as.data.table() %>% 
                    extract(, .(n = sum(count)), by ='index') %>%
                    extract2('n') 
    if (verbose) cmessage( '\t\t\tCount %d-mismatch hits in genome      : %s',
                            mismatch,
                            format(signif(Sys.time() - starttime, 2)))
    matches
}


#' Find target cas9 ranges with no offtargets
#' @param targetranges  GenomicRanges::GRanges
#' @param mismatch max number of mismatches to consider
#' @param verbose  logical(1)
#' @return GenomicRanges::GRanges
#'         mcols(GRanges) contains sequences and match counts:
#'             matches0 = perfect match counts
#'             matches1 = single mismatch counts
#'             matches2 = double mismatch counts
#' @examples
#' \donrun{
#'    require(magrittr)
#'    bedfile <- system.file('extdata/SRF_sites.bed', package = 'multicrispr')
#'    bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#'    targetranges <- read_bed(bedfile, bsgenome)[1:10] %>% flank_fourways()
#'    targetranges %>% find_specific_cas9ranges(mismatch=0)
#' }
#' @export
find_specific_cas9ranges <- function(
    targetranges, 
    mismatch = 2,
    verbose = TRUE
){
    # Assert
    assertive.types::assert_is_a_number(mismatch)

    # Prepare
    bsgenome   <- get_bsgenome(targetranges)
    targetseqs <- seqs(targetranges)
    cas9ranges <- targetranges %>% find_cas9ranges(verbose = verbose)
    cas9seqdt  <- data.table::data.table(
                    seqs = unique(cas9ranges$seqs))
    
    # Count-store-filter for 0-2 mismatches
    if (verbose) message('\tFind cas9seq hits')
    for (mis in 0:mismatch){
        if (verbose) cmessage('\t\twith %d mismatch(es)', mis)
        # Count
        target_matches  <-  cas9seqdt$seqs %>% 
                            count_target_matches(targetseqs, mis, verbose)
        genome_matches  <-  cas9seqdt$seqs %>% 
                            count_genome_matches(bsgenome,   mis, verbose)
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
                        as.granges(bsgenome)
    if (verbose)    cmessage('Returning %d cas9seqs across %d ranges',
                    length(unique(specificranges$seqs)), 
                    length(specificranges))
    specificranges
        
}
