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
#' bedfile  <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' granges <- read_bed(bedfile, bsgenome)
#' sequence(granges)
#' @export
sequence <- function(granges){
    assertive.base::assert_is_identical_to_true(is(granges, 'GRanges'))
    BSgenome::getSeq(get_bsgenome(granges), granges)
}


#=============================================================================
# Find cas9sites
#=============================================================================

#' Find cas9sites in targetranges
#' @param targetranges  GenomicRanges::GRanges
#' @param specific      logical(1)
#' @param verbose       logical(1)
#' @return GenomicRanges::GRanges
#' @examples
#' require(magrittr)
#' bedfile  <- system.file('extdata/SRF_sites.bed', package='crisprapex')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#' targetranges <- read_bed(bedfile, bsgenome) %>% flank_fourways()
#' find_cas9s(targetranges)
#' @export 
find_cas9s <- function(
    targetranges, 
    verbose = TRUE
){
    # Assert
    assertive.base::assert_is_identical_to_true(is(targetranges, 'GRanges'))
    assertive.types::assert_is_a_bool(verbose)
    
    # Find cas9s in targetranges
    if (verbose) message('\tFind N{20}NGG cas9seqs')
    targetdt <- data.table::as.data.table(targetranges)
    targetdt [ , seq := sequence(targetranges) %>% as.character() ]
    res <- targetdt$seq %>% stringi::stri_locate_all_regex('[ACGT]{21}GG')
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
    cas9dt <-   targetdt %>% 
                tidyr::separate_rows(substart, subend) %>%
                data.table::data.table() %>% 
                extract(, substart := as.numeric(substart)) %>% 
                extract(, subend   := as.numeric(subend))
    cas9dt [ , cas9seq   := substr(seq, substart, subend)    ]
    cas9dt [ , seq := NULL ]
    cas9dt [ strand=='+', cas9start := start + substart - 1  ]
    cas9dt [ strand=='+', cas9end   := start + subend   - 1  ]
    cas9dt [ strand=='-', cas9start := end   - subend   + 1  ]
    cas9dt [ strand=='-', cas9end   := end   - substart + 1  ]
    
    cas9ranges  <-  cas9dt [ , list(seqnames = seqnames, 
                                    start    = cas9start, 
                                    end      = cas9end, 
                                    strand   = strand) ] %>% 
                    unique() %>% 
                    as.granges(get_bsgenome(targetranges))
    # Return
    if (verbose)   cmessage('\t\t%d cas9 seqs across %d targetranges', 
                            length(unique(sequence(cas9ranges))), 
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

#' Find target-specific cas9seqs
#' @param targetranges GenomicRanges::GRanges
#' @param verbose      logical(1)
#' @return data.table(cas9seq, 
#'                    matches0, 
#'                    matches1, 
#'                    matches2, 
#'                    seqnames, 
#'                    start, 
#'                    end, 
#'                    width, 
#'                    strand)
#' @examples
#' \donrun{
#'    require(magrittr)
#'    bedfile <- system.file('extdata/SRF_sites.bed', package = 'crisprapex')
#'    bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#'    targetranges <- read_bed(bedfile, bsgenome)[1:10] %>% flank_fourways()
#'    targetranges %>% find_specific_cas9s()
#' }
#' @export
find_specific_cas9s <- function(targetranges, verbose = TRUE){

    # Prepare
    bsgenome   <- get_bsgenome(targetranges)
    targetseqs <- sequence(targetranges)
    cas9ranges <- targetranges %>% find_cas9s(verbose = verbose)
    S4Vectors::mcols(cas9ranges)$cas9seq <- sequence(cas9ranges)
    cas9seqdt  <- data.table::data.table(
                    cas9seq = unique(as.character(cas9ranges$cas9seq)))
    
    # Count-store-filter for 0-2 mismatches
    if (verbose) message('\tFind cas9seq hits')
    for (mismatch in c(0,1,2)){
        if (verbose) cmessage('\t\twith %d mismatch(es)', mismatch)
        # Count
        target_matches  <-  cas9seqdt [ ,   count_target_matches(
                                                cas9seq, 
                                                targetseqs, 
                                                mismatch, 
                                                verbose) ]
        genome_matches  <-  cas9seqdt [ ,   count_genome_matches(
                                                cas9seq, 
                                                bsgenome,   
                                                mismatch, 
                                                verbose) ]
        # Store
        cas9seqdt [ , (sprintf('matches%d', mismatch)) := target_matches ]
        
        # Filter
        assertive.base::assert_all_are_false(genome_matches < target_matches)
        idx <- genome_matches == target_matches
        if (verbose){
            cmessage('\t\t\tKeep %d/%d target-specific cas9seqs', 
                     sum(idx), length(idx), mismatch)
        }
        cas9seqdt %<>% extract(idx)
    }
    
    # Return
    cas9seqdt %>% merge(data.table::as.data.table(cas9ranges), by = 'cas9seq')
}
