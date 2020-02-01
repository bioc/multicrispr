#' Index genome
#' 
#' Bowtie index genome
#' @param bsgenome \code{\link[BSgenome]{BSgenome-class}}
#' @param indexedgenome string: index directory
#' @return string: directory with indexed genome
#' @examples 
#' # index_genome(BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
#' # index_genome(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
#' @export
index_genome <- function(
    bsgenome, 
    indexedgenome = file.path('~/.multicrispr/bowtie', bsgenome@pkgname)
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
#' @return
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
    targets %<>% reduce()
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


run_bowtie <- function(crisprfa, indexdir, outfile, norc){
    
    assertive::assert_all_are_existing_files(crisprfa)
    assertive::assert_all_are_dirs(indexdir)
    assertive::assert_is_a_bool(norc)
    
    Rbowtie::bowtie(
        sequences = crisprfa,
        index     = file.path(indexdir, basename(indexdir)),
        f         = TRUE,  # fasta input
        #m         = 10000, # ignore seqs with m+ alignments
        a         = TRUE,  # report ALL alignments
        v         = 2,     # up to 3 mismatches
        norc      = norc,  # no reverse complement
        outfile   = outfile,
        force     = TRUE)
}

read_bowtie_results <- function(outfile){
    
    assertive::assert_all_are_existing_files(outfile)
    
    dt <- data.table::fread(
            outfile,
            col.names = c('readname', 'strand', 'target', 'position', 
                          'readseq', 'quality', 'matches', 'mismatches'))
    dt[ , mismatch := stringi::stri_count_fixed(mismatches, '>')]
    dt[ ,   list(   MM0 = sum(mismatch==0),
                    MM1 = sum(mismatch==1),
                    MM2 = sum(mismatch==2)),
            by = 'readname' ]
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
#' @param outfile   string: file where to output bowtie results
#' @param count_vars character(3): how to name count variables
#' @param verbose   TRUE (default) or FALSE
#' @seealso \code{\link{index_genome}}, \code{\link{index_targets}}
#' @examples
#' # TFBS example
#' #-------------
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' targets <- extend(bed_to_granges(bedfile, genome = 'mm10'))
#' indexedtargets <- index_targets(targets, bsgenome, indexedtargets)
#' spacers <- find_spacers(targets, bsgenome)
#' crisprseqs <- unique(paste0(spacers$crisprspacer, spacers$crisprpam))
#' match_seqs(crisprseqs, indexedtargets, norc=FALSE)
#' @export
match_seqs <- function(
    seqs, indexdir, norc,
    outfile = paste0(default_outdir(), '/match_', basename(indexdir), '.txt'),
    count_vars = c('MM0', 'MM1', 'MM2'),
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
    run_bowtie(readfa, indexdir, outfile, norc = norc)

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
    readdt[, readname := NULL]
    setnames(readdt, c('MM0', 'MM1', 'MM2'), count_vars)
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
#' @param outfile   string: file where to output bowtie results
#' @param pam       string (default 'NGG') pam pattern to expand
#' @param count_vars character(3): how to name count variables
#' @param verbose   TRUE (default) or FALSE
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
#' @export
match_spacers <- function(
    spacers, 
    indexdir, 
    norc, 
    outfile = paste0(default_outdir(), '/match_', basename(indexdir), '.txt'),
    pam = 'NGG', 
    count_vars = c('MM0', 'MM1', 'MM2'), 
    verbose = TRUE
){
    
    # Expand pams
    if (verbose) cmessage('\t\tExpand iupac ambiguities in pam')
    spacerseqs <- unique(spacers$crisprspacer)
    pamseqs <- expand_iupac_ambiguities(pam)
    crisprdt <- data.table(crisprspacer = rep(spacerseqs, each = length(pamseqs)), 
                        crisprpam     = rep(pamseqs, times = length(spacerseqs))) %>% 
              extract(, crispr := paste0(crisprspacer, pamseqs))
    crisprseqs <- unique(crisprdt$crispr)
    
    # Count matches
    matches <- match_seqs(crisprseqs, 
                            indexdir, 
                            norc       = norc, 
                            outfile    = outfile, 
                            count_vars = c('MM0', 'MM1', 'MM2'), 
                            verbose    = verbose)
    matches %<>% merge(crisprdt, by.x = 'readseq', by.y = 'crispr', all = TRUE)
    matches %<>% extract(, list(MM0 = sum(MM0), MM1 = sum(MM1), MM2 = sum(MM2)), 
                               by = 'crisprspacer')
    matches %>% setnames(c('MM0', 'MM1', 'MM2'), count_vars)
    matches[spacerseqs, on = 'crisprspacer']
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
#'  spacers %>% add_target_counts(targets, bsgenome)
#' @export
add_target_counts <- function(
    spacers, 
    targets,
    bsgenome,
    outdir = default_outdir(),
    pam       = 'NGG',
    verbose   = TRUE
){
    
    # Index targets
    indexedtargets <- paste0(outdir, '/target')
    dir.create(indexedtargets, showWarnings = FALSE, recursive = FALSE)
    indexedtargets <- index_targets(targets, bsgenome, indexedtargets)

    # Match spacers to targets
    if (verbose) cmessage('\tAdd target match counts')
    outfile <- paste0(outdir, '/match_', basename(indexedtargets), '.txt')
    count_vars <- c('T0', 'T1', 'T2')
    matches <- match_spacers(
                    spacers, indexedtargets, norc = TRUE, outfile = outfile,
                    pam = pam, count_vars = count_vars, verbose = verbose)
    
    # Add counts to spacers. Return.
    dt <- as.data.table(spacers) %>% 
        merge(matches, by = 'crisprspacer', sort = FALSE)
    targetvars <- names(dt) %>% extract(stri_startswith_fixed(., 'target'))
    crisprvars <- c('crisprname', 'crisprspacer', 'crisprpam', 'crisprext') %>% 
                    intersect(names(dt))
    othervars <- setdiff(names(dt), c(targetvars, crisprvars))
    dt %<>% extract(, c(targetvars, crisprvars, othervars), with = FALSE)
    GRanges(dt, seqinfo = seqinfo(spacers))
    
}
#' Add genome counts
#' 
#' Count spacer matches to genome and add to GRanges
#' 
#' Expands iupac amgiguities in the pam sequence.
#' Matches all resulting sequences against (indexes) target and genome.
#' Adds match counts to GRanges object, and then returns it.
#' @param spacers      spacer \code{\link[GenomicRanges]{GRanges-class}}
#' @param indexedgenome  dir with bowtie-indexed genome
#' @param outdir       dir where output is written to
#' @param pam          string (default 'NGG') pam pattern to expand
#' @param verbose      TRUE (default) or FALSE
#' @seealso \code{\link{index_genome}}, \code{\link{index_targets}}
#' @examples
#' # PE example
#' #-----------
#'  bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#'  gr <- char_to_granges(c(PRNP = 'chr20:4699600:+',             # snp
#'                          HBB  = 'chr11:5227002:-',             # snp
#'                          HEXA = 'chr15:72346580-72346583:-',   # del
#'                          CFTR = 'chr7:117559593-117559595:+'), # ins
#'                        bsgenome)
#'  spacers <- find_pe_spacers(gr, bsgenome)
#'  indexedgenome <- '~/.multicrispr/bowtie/BSgenome.Hsapiens.UCSC.hg38'
#'  # index_genome(bsgenome, indexedgenome)
#'  # add_genome_counts(spacers, indexedgenome)
#'     
#' # TFBS example
#' #-------------
#'  bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'  bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'  targets <- extend(bed_to_granges(bedfile, genome = 'mm10'))
#'  spacers <- find_spacers(targets, bsgenome)
#'  indexedgenome <- '~/.multicrispr/bowtie/BSgenome.Mmusculus.UCSC.mm10'
#'  # index_genome(bsgenome, indexedgenome)
#'  # add_genome_counts(spacers, genome)
#' @export
add_genome_counts <- function(
    spacers, 
    indexedgenome,
    outdir  = default_outdir(),
    pam     = 'NGG', 
    verbose = TRUE
){
    # Add genome matches
    if (verbose) message('\tAdd genome match counts')
    outfile <- paste0(outdir, '/match_', basename(indexedgenome), '.txt')
    count_vars <- c('G0', 'G1', 'G2')
    matches <- match_spacers(
                        spacers, indexedgenome, norc = FALSE, outfile = outfile, 
                        pam = pam, count_vars = count_vars, verbose = verbose)
    dt <- spacers %>% 
        as.data.table() %>% 
        merge(matches, by = 'crisprspacer', sort = FALSE)
    
    # Organize columns
    targetvars <- names(dt) %>% extract(stringi::stri_startswith_fixed(., 'target'))
    crisprvars <- c('crisprname', 'crisprspacer', 'crisprpam', 'crisprext') %>% 
                    intersect(names(dt))
    othervars <- setdiff(names(dt), c(targetvars, crisprvars))
    dt %<>% extract(, c(targetvars, crisprvars, othervars), with = FALSE)
    
    # Return GRanges
    GRanges(dt, seqinfo = seqinfo(spacers))
}

default_outdir <- function() paste0(tempdir(), '/multicrispr/bowtie')



#' Filter for target-specific spacers
#' 
#' @param spacers   spacer \code{\link[GenomicRanges]{GRanges-class}}
#' @param targets   target \code{\link[GenomicRanges]{GRanges-class}}
#' @param indexedgenome  bowtie indexed genome dir
#' @param outdir    directory where output is written to
#' @param pam       string (default 'NGG'): pam sequence
#' @param verbose   TRUE (default) or FALSE
#' @examples
#' # TFBS example
#' #-------------
#'  bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'  bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'  targets <- extend(bed_to_granges(bedfile, genome = 'mm10'))
#'  spacers <- find_spacers(targets, bsgenome)
#'  indexedgenome <- '~/.multicrispr/bowtie/BSgenome.Mmusculus.UCSC.mm10'
#'  # index_genome(indexedgenome)
#'  # filter_target_specific(spacers, targets, bsgenome, indexedgenome)
#' @export
filter_target_specific <- function(
    spacers,
    targets,
    bsgenome, 
    indexedgenome, 
    outdir = default_outdir(), 
    pam       = 'NGG', 
    verbose   = TRUE
){
    # Add target matches
    spacers %<>% add_target_counts(
                    targets, bsgenome, outdir, pam=pam, verbose=verbose)
    
    # Add genome matches
    spacers %<>% add_genome_counts(indexedgenome, outdir, pam=pam, verbose=verbose)
    
    # Sanity check
    assertive::assert_is_of_length(subset(spacers, T0>G0), 0)
    assertive::assert_is_of_length(subset(spacers, T1>G1), 0)
    assertive::assert_is_of_length(subset(spacers, T2>G2), 0)
    
    # Filter for specificity
    digits <- ceiling(log10(length(spacers)))
    if (verbose) cmessage('\tFilter %d spacers', length(spacers))
    for (mis in 0:2){
        idx <-  mcols(spacers)[[paste0('T', mis)]] == 
                mcols(spacers)[[paste0('G', mis)]]
        spacers %<>% subset(idx)
        if (verbose) cmessage('\t       %s T%d==G%d', 
                                format(sum(idx), width = digits), mis, mis)
    }
    
    # Return
    spacers
}


#' Filter for PE-specific spacers
#' 
#' Filters for spacers which hit a single prime editing site
#' 
#' @param spacers   spacer \code{\link[GenomicRanges]{GRanges-class}}
#' @param targets   target \code{\link[GenomicRanges]{GRanges-class}}
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
    spacers, genomedir, outdir = default_outdir(), pam = 'NGG', verbose = TRUE
){
    # Add genome matches
    spacers %<>% add_genome_counts(genomedir, outdir, pam=pam, verbose=verbose)
    
    # Filter for specificity
    digits <- ceiling(log10(length(spacers)))
    if (verbose) cmessage('\tFilter %d spacers', length(spacers))
    idx <- spacers$G0==1
    spacers %<>% subset(idx)
    if (verbose) cmessage('\t       %s G0==1', format(sum(idx), width = digits))
    
    # Return
    spacers
}

