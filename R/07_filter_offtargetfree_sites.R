#============================================================================

#' Count target/genome matches
#' @param crisprseqs  character vector or \code{\link[Biostrings]{XStringSet-class}}
#' @param targetseqs  character vector or \code{\link[Biostrings]{XStringSet-class}}
#' @param bsgenome    \code{\link[BSgenome]{BSgenome-class}}
#' @param mismatch    number: number of allowed mismatches 
#' @param include     character vector: offtarget analysis chromosomes
#' @param verbose     TRUE (default) or FALSE
#' @return numeric(length(crisprseqs))
#' @examples
#' # Read ranges and extend
#'     require(magrittr)
#'     bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'     targets <- bed_to_granges(bedfile, 'mm10')
#'     targets %<>% extend()
#' 
#' # Add seqs and find crispr sites
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'     targets %<>% add_seq(bsgenome)
#'     sites <- find_crispr_sites(targets)
#' 
#' # Count matches
#'     crisprseqs <- sites$seq[1:10]
#'     count_target_matches(crisprseqs, targets$seq, 0)
#'     count_genome_matches(crisprseqs, bsgenome,    0, include = 'chrY')
#' @export
count_target_matches <- function(
    crisprseqs, 
    targetseqs, 
    mismatch, 
    verbose = TRUE
){
    
    # Assert
    assertive.types::assert_is_any_of(crisprseqs, c('character', 'XStringSet'))
    assertive.types::assert_is_any_of(targetseqs, c('character', 'XStringSet'))
    assertive.types::assert_is_a_number(mismatch)
    assertive.sets::assert_is_subset(mismatch, c(0,1,2))
    assertive.types::assert_is_a_bool(verbose)
    
    # Count
    uniquecas9s <- unique(crisprseqs)
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
    unname(matches[crisprseqs])
}


#' Convert "included" into "excluded" seqnames
#'
#' @details Convert "included" into "excluded" seqnames. 
#'          Do this in regex format as required by 
#' \code{\link[BSgenome]{BSgenome-utils}}) \code{vcountPDict} and 
#' \code{vmatchPDict}
#' @param bsgenome \code{\link[BSgenome]{BSgenome-class}}
#' @param include character vector with names of included chromosomes
#' @examples 
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#' include_to_exclude(bsgenome)
#' @export
include_to_exclude <- function(
    bsgenome, 
    include = GenomeInfoDb::standardChromosomes(bsgenome)
){
    assertive.types::assert_is_all_of(bsgenome, 'BSgenome')
    assertive.types::assert_is_character(include)
    
    # needs to be a regex: https://support.bioconductor.org/p/126593/
    sprintf('^%s$', setdiff(BSgenome::seqnames(bsgenome), include))
}


#' @rdname count_target_matches
#' @export
count_genome_matches <- function(
    crisprseqs, 
    bsgenome, 
    mismatch,
    include = BSgenome::seqnames(bsgenome),
    verbose     = TRUE
){
    # Assert
    assertive.types::assert_is_any_of(crisprseqs, c('character', 'XStringSet'))
    assertive.types::assert_is_any_of(bsgenome, 'BSgenome')
    assertive.sets::assert_is_subset(mismatch, c(0,1,2))
    assertive.types::assert_is_character(include)
    assertive.types::assert_is_a_bool(verbose)

    # Comply
    . <- count <- NULL

    # Count
    starttime <- Sys.time()
    uniquecas9s <- unique(crisprseqs)
    matches  <- Biostrings::vcountPDict(
                    Biostrings::DNAStringSet(uniquecas9s),
                    bsgenome,
                    min.mismatch = mismatch,
                    max.mismatch = mismatch, 
                    exclude      = include_to_exclude(bsgenome, include), 
                    verbose      = verbose) %>% 
                data.table::as.data.table() %>% 
                magrittr::extract(, .(n = sum(count)), by ='index') %>%
                magrittr::extract2('n') %>% 
                magrittr::set_names(uniquecas9s)
    
    # Return
    if (verbose) cmessage( '\t\t\tCount %d-mismatch hits in genome      : %s',
                            mismatch, format(signif(Sys.time() - starttime, 2)))
    unname(matches[crisprseqs])
}


add_seqinfo <- function(gr, bsgenome){
    GenomeInfoDb::seqinfo(gr) <- GenomeInfoDb::seqinfo(bsgenome)
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
#'     targets %<>% extract(GenomicRanges::seqnames(.)=='chrY')
#' 
#' # Add seqs and find crispr sites
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'     targets %<>% add_seq(bsgenome)
#'     sites <- find_crispr_sites(targets)
#'    
#' # Filter for offtarget-free sites
#'     filter_offtargetfree_sites(
#'         sites, targets, bsgenome, 0, include = 'chrY')
#' @export
filter_offtargetfree_sites <- function(sites, targets, bsgenome, mismatch = 2, 
    include = GenomeInfoDb::standardChromosomes(targets), 
    plot = TRUE, verbose = TRUE){
    
    # Assert. Prepare
    assertive.types::assert_is_all_of(sites,    'GRanges')
    assertive.types::assert_is_all_of(targets,  'GRanges')
    assertive.types::assert_is_all_of(bsgenome, 'BSgenome')
    assertive.types::assert_is_a_number(mismatch)
    assertive.types::assert_is_character(include)
    assertive.sets::assert_is_subset(include, BSgenome::seqnames(bsgenome))
    assertive.sets::assert_is_subset('seq',names(GenomicRanges::mcols(sites)))
    assertive.sets::assert_is_subset('seq',names(GenomicRanges::mcols(targets)))
    targetseqs <- targets$seq
    cas9seqdt  <- data.table::data.table(seq = unique(sites$seq))
    
    # Count-store-filter for 0-2 mismatches
    if (verbose) message('\tFind crispr seq (mis)matches')
    for (mis in 0:mismatch){
        if (verbose) cmessage('\t\twith %d mismatch(es)', mis)
        # Count
        target_matches <- count_target_matches(
                            cas9seqdt$seq, targetseqs, mis, verbose = verbose)
        genome_matches <- count_genome_matches(cas9seqdt$seq, bsgenome, mis, 
                            include = include, verbose = verbose)
        # Store
        cas9seqdt [ , (sprintf('matches%d', mis)) := target_matches ]
        
        # Filter
        assertive.base::assert_all_are_false(genome_matches < target_matches)
        idx <- genome_matches == target_matches
        if (verbose)   cmessage('\t\t\tKeep %d/%d offtarget-free crispr seqs', 
                                sum(idx), length(idx), mis)
        cas9seqdt %<>% magrittr::extract(idx)
    }
    
    # Filter for offtargetfree sites
    offtargetfree  <-   cas9seqdt %>% 
                        merge(data.table::as.data.table(sites), by = 'seq') %>% 
                        methods::as('GRanges') %>%  add_seqinfo(bsgenome)
    # Plot/Report
    if (plot)    plot_karyogram(GenomicRanges::GRangesList(
                                    target = targets, cas9 = sites, 
                                    offtargetfree = offtargetfree))
    if (verbose) cmessage('\tReturn %d crispr seqs across %d ranges',
                    length(unique(offtargetfree$seq)), length(offtargetfree))
    # Return
    offtargetfree
}
