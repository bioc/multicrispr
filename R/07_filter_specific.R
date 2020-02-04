default_outdir <- function() paste0(tempdir(), '/multicrispr/bowtie')

default_indexedgenome <- function(x){
    bsgenome <- if (methods::is(x, 'BSgenome')){ 
                    x 
    } else if (methods::is(x, 'GRanges')){ 
                    getBSgenome(genome(x)[1])}
    file.path('~/.multicrispr/bowtie', bsgenome@pkgname)
}

#' Index genome
#' 
#' Bowtie index genome
#' @param bsgenome \code{\link[BSgenome]{BSgenome-class}}
#' @param indexedgenome string: directory with bowtie-indexed genome
#' @return string: directory with indexed genome
#' @examples 
#' # index_genome(BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
#' # index_genome(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
#' @export
index_genome <- function(
    bsgenome, 
    indexedgenome = default_indexedgenome(bsgenome)
){

    # Assert
    assert_is_all_of(bsgenome, 'BSgenome')

    # Create Names
    dir.create(indexedgenome, showWarnings = FALSE, recursive = TRUE)
    genomefa <- paste0(dirname(indexedgenome), '/', bsgenome@pkgname, '.fa')

    # Write to fasta
    BSgenome::writeBSgenomeToFasta(bsgenome, genomefa)

    # Index genome
    subdirs <- list.dirs(indexedgenome, full.names = FALSE, recursive = FALSE)
    Rbowtie::bowtie_build(  genomefa,
                            indexedgenome,
                            prefix = basename(indexedgenome), 
                            force = TRUE)
    
    # Return
    return(indexedgenome)
}


#' Index targets
#' 
#' Bowtie index targets
#' @param targets   \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome  \code{\link[BSgenome]{BSgenome-class}}
#' @param indexedtargets    string: output directory
#' @param verbose           TRUE (default) or FALSE 
#' @return indextargets: output directory
#' @examples 
#' require(magrittr)
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' targets <- extend(bed_to_granges(bedfile, genome = 'mm10'))
#' index_targets(targets, bsgenome)
#' @export
index_targets <- function(
    targets, 
    bsgenome,
    indexedtargets = paste0(default_outdir(), '/targets'), 
    verbose = TRUE
){
    # Assert
    assert_is_all_of(targets,  'GRanges')

    # Add inverse strand
    if (verbose) cmessage('\tIndex target sequences')
    if (verbose) cmessage('\t\t%d target ranges', length(targets))
    targets %<>% add_inverse_strand(plot = FALSE, verbose = TRUE)
    targets %<>% sort(ignore.strand = TRUE)

    # Reduce
    targets %<>% GenomicRanges::reduce()
    if (verbose) cmessage('\t\t%d ranges after merging overlaps', length(targets))
    
    # Write to fasta
    targetfa <- file.path(dirname(indexedtargets), 'target.fa')
    targetseqs <- BSgenome::getSeq(bsgenome, targets)
    names(targetseqs) <- sprintf(
                          '%s:%s-%s:%s',
                          as.character(seqnames(targets)),
                          start(targets),
                          end(targets),
                          strand(targets))
    if (verbose) cmessage('\t\tWrite seqs  to %s', targetfa)
    dir.create(dirname(targetfa), showWarnings = FALSE, recursive = TRUE)
    Biostrings::writeXStringSet(targetseqs, targetfa)

    # Index targets
    if (verbose) cmessage('\t\tWrite index to %s', indexedtargets)
    Rbowtie::bowtie_build(  targetfa,
                            indexedtargets,
                            prefix = basename(indexedtargets),
                            force = TRUE)
    
    # Return
    return(indexedtargets)
}


run_bowtie <- function(crisprfa, indexdir, outfile, norc, mismatches = 2){
    
    assertive::assert_all_are_existing_files(crisprfa)
    assertive::assert_all_are_dirs(indexdir)
    assertive::assert_is_a_bool(norc)
    assertive::assert_is_a_number(mismatches)
    assertive::assert_is_subset(mismatches, 0:3)
    
    Rbowtie::bowtie(
        sequences = crisprfa,
        index     = file.path(indexdir, basename(indexdir)),
        f         = TRUE,        # fasta input
        #m         = 10000,      # ignore seqs with m+ alignments
        a         = TRUE,        # report ALL alignments
        v         = mismatches,  # up to 3 mismatches
        norc      = norc,        # no reverse complement
        outfile   = outfile,
        force     = TRUE)
}

read_bowtie_results <- function(outfile){
    
    assertive::assert_all_are_existing_files(outfile)
    mismatch <- mismatches <- NULL
    
    dt <- data.table::fread(
            outfile,
            col.names = c('readname', 'strand', 'target', 'position', 
                          'readseq', 'quality', 'matches', 'mismatches'))
    dt[ , mismatch := stringi::stri_count_fixed(mismatches, '>')]
    counts <- data.table::data.table(readname = unique(dt$readname))
    for (mis in sort(unique(dt$mismatch))){
        counts %<>% merge( dt[ , list(counts = sum(mismatch==mis)), by = 'readname'] )
        counts %>% data.table::setnames('counts', paste0('MM', mis))
    }
    counts
}


#' Add match counts
#' 
#' Count matches to indexed target/genome and add to GRanges
#' 
#' \code{match_seqs} matches sequences against (indexed) target and genome
#' \code{add_crispr_matches} expands iupac amgiguities in the pam sequence, 
#' matches all resulting sequences against (indexes) target and genome, 
#' adds match counts to GRanges object, and then returns it
#' 
#' @param seqs      character vector: sequences to match against indexed ref
#' @param indexdir  string: dir containing indexed reference.
#'                  This can be an indexed genome( \code{\link{index_genome}}
#'                  It can also be indexed targets (\code{\link{index_targets}})
#' @param norc      TRUE or FALSE: whether to run bowtie also with revcompls
#'                  Generally TRUE for genome and FALSE for target matches, 
#'                  because target ranges generally include both strands.
#' @param mismatches max number of mismatches to consider
#' @param outfile   string: file where to output bowtie results
#' @param verbose   TRUE (default) or FALSE
#' @seealso \code{\link{index_genome}}, \code{\link{index_targets}}
#' @examples
#' # TFBS example
#' #-------------
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' targets <- extend(bed_to_granges(bedfile, genome = 'mm10'))
#' indexedtargets <- index_targets(targets, bsgenome)
#' spacers <- find_spacers(targets, bsgenome)
#' crisprseqs <- unique(paste0(spacers$crisprspacer, spacers$crisprpam))
#' match_seqs(crisprseqs, indexedtargets, norc=FALSE)
#' match_seqs(crisprseqs, indexedtargets, norc=FALSE, mismatch=3)
#' @export
match_seqs <- function(
    seqs, indexdir, norc, mismatches = 2,
    outfile = paste0(default_outdir(), '/match_', basename(indexdir), '.txt'),
    verbose    = TRUE
){

    # Assert
    assertive::assert_is_character(seqs)
    assertive::assert_has_no_duplicates(seqs)

    # Write reads to fasta
    reads <- Biostrings::DNAStringSet(unique(seqs))
    reads %<>% name_uniquely('read')
    readfa <- file.path(dirname(outfile), 'reads.fa')
    dir.create(dirname(readfa), recursive = TRUE, showWarnings = FALSE)
    if (verbose) cmessage('\t\tWrite reads to %s', readfa)
    Biostrings::writeXStringSet(reads, readfa)

    # Map reads
    if (verbose) cmessage('\t\tMap reads: %s', outfile)
    run_bowtie(readfa, indexdir, outfile, norc = norc, mismatches = mismatches)

    # Load results
    if (verbose) cmessage('\t\tLoad results')
    matches <- read_bowtie_results(outfile)
    
    # Count matches
    if (verbose) cmessage('\t\tCount matches')
    readdt <- data.table(readname = names(reads), 
                        readseq   = unname(as.character(reads)))
    readdt %<>% merge(matches, by='readname', all=TRUE, sort=FALSE)

    # Turn NA into 0
    readdt <- cbind(readdt[, 1:2], 
                    setnafill(readdt[, 3:ncol(readdt)], fill = 0))

    # Return
    readdt[, 'readname' := NULL]
    readdt[seqs, on = 'readseq']
}

explode <- function(x) unlist(strsplit(x, character(0)))
paste_dtcols <- function(x) x [, do.call(paste0, .SD) ]
expand_iupac_ambiguities <- function(x){
    paste_dtcols( 
        as.data.table( 
            lapply(Biostrings::IUPAC_CODE_MAP[explode(x)], explode)))
}

#' Match spacers
#' 
#' Count matches to indexed target/genome and add to GRanges
#' 
#' Expands iupac amgiguities in the pam sequence.
#' Matches all resulting sequences against (indexes) target and genome.
#' Adds match counts to GRanges object, and then returns it.
#' @param spacers   spacer \code{\link[GenomicRanges]{GRanges-class}}
#' @param indexdir  string: dir containing indexed reference.
#'                  This can be an indexed genome( \code{\link{index_genome}}
#'                  It can also be indexed targets (\code{\link{index_targets}})
#' @param norc      TRUE or FALSE: whether to run bowtie also with revcompls
#'                  Generally TRUE for genome and FALSE for target matches, 
#'                  because target ranges generally include both strands.
#' @param mismatches number (default 2): max number of mismatches to consider
#' @param outfile    string: file where to output bowtie results
#' @param pam        string (default 'NGG') pam pattern to expand
#' @param verbose    TRUE (default) or FALSE
#' @seealso \code{\link{index_genome}}, \code{\link{index_targets}}
#' @examples
#' # TFBS example
#' #-------------
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' targets <- extend(bed_to_granges(bedfile, genome = 'mm10'))
#' indexedtargets <- index_targets(targets, bsgenome, indexedtargets)
#' spacers <- find_spacers(targets, bsgenome)
#' match_spacers(spacers, indexedtargets, norc=FALSE)
#' match_spacers(spacers, indexedtargets, norc=FALSE, mismatches = 3)
#' @export
match_spacers <- function(
    spacers, 
    indexdir, 
    norc, 
    mismatches = 2,
    outfile = paste0(default_outdir(), '/match_', basename(indexdir), '.txt'),
    pam     = 'NGG', 
    verbose = TRUE
){
    # Comply
    crispr <- crisprspacer <- . <- NULL
    
    # Expand pams
    if (verbose) cmessage('\t\tExpand iupac ambiguities in pam')
    spacerseqs <- unique(spacers$crisprspacer)
    pamseqs <- expand_iupac_ambiguities(pam)
    crisprdt <- data.table(
                    crisprspacer = rep(spacerseqs, each = length(pamseqs)), 
                    crisprpam    = rep(pamseqs, times = length(spacerseqs))) %>% 
              extract(, crispr := paste0(crisprspacer, pamseqs))
    crisprseqs <- unique(crisprdt$crispr)
    
    # Count matches
    matches <- match_seqs(crisprseqs, 
                            indexdir, 
                            norc       = norc, 
                            mismatches = mismatches,
                            outfile    = outfile, 
                            verbose    = verbose)
    matches %<>% merge(crisprdt, by.x = 'readseq', by.y = 'crispr', all = TRUE)
    
    # Aggregate per spacer
    count_vars <- names(matches) %>% extract(stri_startswith_fixed(., 'MM'))
    spacer_dt <- data.table(crisprspacer = unique(matches$crisprspacer))
    for (var in count_vars){
        spacer_dt %<>% merge(matches[ , list(counts = sum(get(var))), 
                                     by = 'crisprspacer'], 
                            by = 'crisprspacer')
        setnames(spacer_dt, 'counts', var)
    }
    
    # Return
    spacer_dt[spacerseqs, on = 'crisprspacer']
}

#' Add target counts
#' 
#' Count spacer matches among targets
#' 
#' Expands iupac amgiguities in the pam sequence.
#' Matches all resulting sequences against (indexes) target and genome.
#' Adds match counts to GRanges object, and then returns it.
#' @param spacers  spacer \code{\link[GenomicRanges]{GRanges-class}}
#' @param targets  target \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome \code{\link[BSgenome]{BSgenome-class}}
#' @param mismatches number (default 2): max number of mismatches to consider
#' @param outdir   dir where output is written to
#' @param pam      string (default 'NGG') pam pattern to expand
#' @param verbose  TRUE (default) or FALSE
#' @seealso \code{\link{index_genome}}, \code{\link{index_targets}}
#' @examples
#' # TFBS example
#' #-------------
#'  bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'  genome <- '~/.multicrispr/bowtie/BSgenome.Mmusculus.UCSC.mm10'
#'  bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'  targets  <- extend(bed_to_granges(bedfile, genome = 'mm10'))
#'  spacers <- find_spacers(targets, bsgenome)
#'  add_target_counts(spacers, targets, bsgenome)
#' @export
add_target_counts <- function(
    spacers, 
    targets,
    bsgenome,
    mismatches = 2,
    outdir     = default_outdir(),
    pam        = 'NGG',
    verbose    = TRUE
){
    # Comply
    . <- NULL
    
    # Index targets
    indexedtargets <- paste0(outdir, '/target')
    dir.create(indexedtargets, showWarnings = FALSE, recursive = FALSE)
    indexedtargets <- index_targets(targets, bsgenome, indexedtargets)

    # Match spacers to targets
    if (verbose) cmessage('\tAdd target match counts')
    outfile <- paste0(outdir, '/match_', basename(indexedtargets), '.txt')
    matches <- match_spacers(
                    spacers, indexedtargets, norc = TRUE, 
                    mismatches = mismatches, outfile = outfile,
                    pam = pam, verbose = verbose)
    names(matches) %<>% stringi::stri_replace_first_regex('^MM', 'T')
    
    # Add counts to spacers
    dt  <-  gr2dt(spacers) %>%
            merge(matches, by = 'crisprspacer', sort = FALSE)
    
    # Organize columns and return
    targetvars <- names(dt) %>% extract(stri_startswith_fixed(., 'target'))
    crisprvars <- c('crisprname', 'crisprspacer', 'crisprpam', 'crisprext') %>% 
                    intersect(names(dt))
    othervars <- setdiff(names(dt), c(targetvars, crisprvars))
    
    dt %>% 
    extract(, c(targetvars, crisprvars, othervars), with = FALSE) %>% 
    dt2gr(seqinfo(spacers))
    
}
#' Add genome counts
#' 
#' Count spacer matches to genome and add to GRanges
#' 
#' Expands iupac amgiguities in the pam sequence.
#' Matches all resulting sequences against (indexes) target and genome.
#' Adds match counts to GRanges object, and then returns it.
#' @param spacers        spacer \code{\link[GenomicRanges]{GRanges-class}}
#' @param indexedgenome  dir with bowtie-indexed genome
#' @param mismatches  number (default 2): max number of mismatches to consider
#' @param outdir         dir where output is written to
#' @param pam            string (default 'NGG') pam pattern to expand
#' @param verbose        TRUE (default) or FALSE
#' @seealso \code{\link{index_genome}}, \code{\link{index_targets}}
#' @examples
#' # PE example
#' #-----------
#'  require(magrittr)
#'  bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#'  gr <- char_to_granges(c(PRNP = 'chr20:4699600:+',             # snp
#'                          HBB  = 'chr11:5227002:-',             # snp
#'                          HEXA = 'chr15:72346580-72346583:-',   # del
#'                          CFTR = 'chr7:117559593-117559595:+'), # ins
#'                        bsgenome)
#'  spacers <- find_pe_spacers(gr, bsgenome)
#'  indexedgenome <- '~/.multicrispr/bowtie/BSgenome.Hsapiens.UCSC.hg38'
#'  # index_genome(bsgenome, indexedgenome)
#'  # add_genome_counts(spacers, indexedgenome, mismatches=1)
#'     
#' # TFBS example
#' #-------------
#'  bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'  bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'  targets <- extend(bed_to_granges(bedfile, genome = 'mm10'))
#'  spacers <- find_spacers(targets, bsgenome)
#'  indexedgenome <- '~/.multicrispr/bowtie/BSgenome.Mmusculus.UCSC.mm10'
#'  # index_genome(bsgenome, indexedgenome)
#'  # add_genome_counts(spacers, indexedgenome)
#'  # add_genome_counts(spacers, indexedgenome, mismatches=3)
#' @export
add_genome_counts <- function(
    spacers, 
    indexedgenome = default_indexedgenome(spacers),
    mismatches    = 2,
    outdir        = default_outdir(),
    pam           = 'NGG', 
    verbose       = TRUE
){
    . <- NULL
    
    # Add genome matches
    if (verbose) message('\tAdd genome match counts')
    outfile <- paste0(outdir, '/match_', basename(indexedgenome), '.txt')
    matches <- match_spacers(
                        spacers, indexedgenome, norc = FALSE, 
                        mismatches = mismatches, outfile = outfile, 
                        pam = pam, verbose = verbose)
    names(matches) %<>% stringi::stri_replace_first_regex('^MM', 'G')
    dt <- gr2dt(spacers) %>%
            merge(matches, by = 'crisprspacer', sort = FALSE)
    
    # Organize columns and return
    targetvars <- names(dt) %>% extract(stringi::stri_startswith_fixed(., 'target'))
    crisprvars <- c('crisprname', 'crisprspacer', 'crisprpam', 'crisprext') %>% 
                    intersect(names(dt))
    othervars <- setdiff(names(dt), c(targetvars, crisprvars))
    
    dt %>% 
    extract(, c(targetvars, crisprvars, othervars), with = FALSE) %>% 
    dt2gr(seqinfo(spacers))
}


#' @rdname filter_target_specific
#' @export
add_specificity <- function(
    spacers,
    targets,
    bsgenome, 
    indexedgenome = default_indexedgenome(spacers), 
    mismatches    = 2,
    outdir        = default_outdir(), 
    pam           = 'NGG',
    plot          = TRUE,
    verbose       = TRUE
){
    # Comply
    specific <- NULL
    
    # Add target/genome matches
    spacers %<>% add_target_counts(
                    targets, bsgenome, mismatches = mismatches, outdir = outdir,
                    pam = pam, verbose = verbose)
    spacers %<>% add_genome_counts(
                    indexedgenome, mismatches = mismatches, outdir = outdir, 
                    pam = pam, verbose = verbose)
    
    # Sanity check
    for (mis in 0:mismatches){
            assertive.base::assert_all_are_true(
                mcols(spacers)[[paste0('T', mis)]] <=
                mcols(spacers)[[paste0('G', mis)]])
    }

    # Add specificity info
    digits <- ceiling(log10(length(spacers)))
    if (verbose) cmessage('\tFilter %d spacers', length(spacers))
    spacers$specific <- FALSE
    idx <- rep(TRUE, length(spacers))
    for (mis in 0:2){
        idx %<>% and(mcols(spacers)[[paste0('T', mis)]] == 
                     mcols(spacers)[[paste0('G', mis)]])
        if (verbose) cmessage('\t       %s T%d==G%d', 
                                format(sum(idx), width = digits), mis, mis)
    }
    spacers$specific <- idx
    
    # Plot
    if (plot){
        grplot <- gr2dt(spacers) %>% 
                    extract(, if(any(specific)) .SD, by = 'targetname') %>% 
                    dt2gr(seqinfo(spacers))
        p <- plot_intervals(grplot, alpha_var = 'specific') + 
            ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.25))
        print(p)
    }
    
    # Return
    spacers
}


#' Filter for target-specific spacers
#' 
#' @param spacers   spacer \code{\link[GenomicRanges]{GRanges-class}}
#' @param targets   target \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome  \code{\link[BSgenome]{BSgenome-class}}
#' @param indexedgenome  bowtie indexed genome dir
#' @param mismatches number (default 2): max number of mismatches to consider
#' @param outdir    directory where output is written to
#' @param pam       string (default 'NGG'): pam sequence
#' @param plot      TRUE (default) or FALSE
#' @param verbose   TRUE (default) or FALSE
#' @examples
#' # TFBS example
#' #-------------
#'  require(magrittr)
#'  bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'  bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'  targets <- extend(bed_to_granges(bedfile, genome = 'mm10'))
#'  spacers <- find_spacers(targets, bsgenome)
#'  indexedgenome <- '~/.multicrispr/bowtie/BSgenome.Mmusculus.UCSC.mm10'
#'  # index_genome(indexedgenome)
#'  # spacers %<>% add_specificity(targets, bsgenome, indexedgenome)
#'  # spacers %>% filter_target_specific(targets, bsgenome, indexedgenome)
#' @export
filter_target_specific <- function(
    spacers,
    targets,
    bsgenome, 
    indexedgenome = default_indexedgenome(spacers), 
    mismatches    = 2,
    outdir        = default_outdir(), 
    pam           = 'NGG',
    plot          = TRUE,
    verbose       = TRUE
){
    # Add specificty info
    spacers %<>% add_specificity(
                    targets       = targets, 
                    bsgenome      = bsgenome, 
                    indexedgenome = indexedgenome, 
                    mismatches    = mismatches,
                    outdir        = outdir, 
                    pam           = pam, 
                    plot          = plot,
                    verbose       = verbose)
    
    # Subset
    specific <- NULL
    spacers %<>% subset(specific==TRUE)
    spacers$specific <- FALSE
    
    # Return
    return(spacers)
}



#' Filter for PE-specific spacers
#' 
#' Filters for spacers which hit a single prime editing site
#' 
#' @param spacers   spacer \code{\link[GenomicRanges]{GRanges-class}}
#' @param indexedgenome  bowtie indexed genome dir
#' @param outdir    directory where output is written to
#' @param pam       string (default 'NGG'): pam sequence
#' @param verbose   TRUE (default) or FALSE
#' @examples
#' # PE example
#' #-----------
#'  require(magrittr)
#'  bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#'  gr <- char_to_granges(c(PRNP = 'chr20:4699600:+',             # snp
#'                          HBB  = 'chr11:5227002:-',             # snp
#'                          HEXA = 'chr15:72346580-72346583:-',   # del
#'                          CFTR = 'chr7:117559593-117559595:+'), # ins
#'                        bsgenome)
#'  spacers <- find_pe_spacers(gr, bsgenome)
#'  indexedgenome <- '~/.multicrispr/bowtie/BSgenome.Hsapiens.UCSC.hg38'
#'  # index_genome(bsgenome, indexedgenome) # one time effort - takes few h
#'  # filter_prime_specific(spacers, indexedgenome)
#' @export
filter_prime_specific <- function(
    spacers, 
    indexedgenome = default_indexedgenome(spacers),
    outdir        = default_outdir(), 
    pam           = 'NGG', 
    verbose       = TRUE
){
    # Add genome matches
    spacers %<>% add_genome_counts(
                    indexedgenome, mismatches = 1, outdir = outdir, 
                    pam = pam, verbose = verbose)
    
    # Filter for specificity
    digits <- ceiling(log10(length(spacers)))
    if (verbose) cmessage('\tFilter %d spacers', length(spacers))
    idx <- spacers$G0==1
    spacers %<>% subset(idx)
    if (verbose) cmessage('\t       %s G0==1', format(sum(idx), width = digits))
    
    # Return
    spacers
}

