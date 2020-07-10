INDEXEDGENOMESDIR <- '~/multicrisprout/indexedgenomes'
    genome_dir <- function(indexedgenomesdir = INDEXEDGENOMESDIR, bsgenome){
            paste0(indexedgenomesdir, '/', bsgenome@pkgname)}
    
    genome_fasta <- function(indexedgenomesdir = INDEXEDGENOMESDIR, bsgenome){ 
            paste0(indexedgenomesdir, '/', bsgenome@pkgname, '.fa')}

OUTDIR <- '~/multicrisprout'
    target_dir      <- function(outdir = OUTDIR){
            paste0(outdir,  '/targets') }
    
    target_fasta    <- function(outdir = OUTDIR){
            paste0(outdir, '/targets.fa')}
    
    spacer_matchfile <- function(outdir = OUTDIR, indexdir){
            paste0(outdir, '/spacers/spacers_to_', basename(indexdir), '.txt')}
    
    spacer_fasta <- function(outdir = OUTDIR){
            paste0(outdir, '/spacers.fa') }

#' Has been indexed?
#' @param bsgenome BSgenome
#' @param indexedgenomesdir directory with indexed genomes
#' @return TRUE or FALSE
#' @examples 
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#' has_been_indexed(bsgenome)
#' @export
has_been_indexed <- function(bsgenome, indexedgenomesdir = INDEXEDGENOMESDIR){
    dir.exists(genome_dir(indexedgenomesdir, bsgenome = bsgenome))
    # don't rm this function - is used in vignette
}


indexed_genomes_s3 <- c("BSgenome.Hsapiens.UCSC.hg38", 
                        "BSgenome.Mmusculus.UCSC.mm10")

#' Index genome
#' 
#' Bowtie index genome
#' 
#' Checks whether already available locally. If not, checks whether indexed
#' version can be downloaded from our s3 storage. If not, builds the
#' index with bowtie. This can take a few hours, but is a one-time operation.
#' @param bsgenome \code{\link[BSgenome]{BSgenome-class}}
#' @param indexedgenomesdir  string: directory with bowtie-indexed genome
#' @param download TRUE (default) or FALSE: whether to download pre-indexed 
#'                 version if available
#' @param overwrite TRUE or FALSE (default)
#' @return invisible(genomdir)
#' @examples
#' bsgenome <- BSgenome.Scerevisiae.UCSC.sacCer1::Scerevisiae
#' index_genome(bsgenome, indexedgenomesdir = tempdir())
#' @export
index_genome <- function(
    bsgenome, indexedgenomesdir = INDEXEDGENOMESDIR, 
    download = TRUE, overwrite = FALSE
){
    # Assert
    assert_is_all_of(bsgenome, 'BSgenome')
    bsname <- GenomeInfoDb::bsgenomeName(bsgenome)

    # Return if already available
    genomedir    <- genome_dir(  indexedgenomesdir, bsgenome)
    message('\tIndex ', bsname)
    if (!overwrite & 
        dir.exists(genomedir) & length(list.files(genomedir))!=0){
        message('\t\tNot required: ', genomedir, ' already exists', '
                Set overwrite=TRUE to overwrite')
        return(invisible(genomedir))
    }

    # Download if present on s3 storage
    if (download & (bsname %in% indexed_genomes_s3)){
        url <- paste0('https://s3.mpi-bn.mpg.de/data-multicrispr-2020/', 
                    bsname, '.zip')
        zipfile <- paste0(genomedir, '.zip')
        message('\t\tDownload pre-indexed version
                For a fresh build instead, set download = FALSE')
        res <- tryCatch(download.file(url, zipfile), error = function(e) 1)
        if (res==0){
            message('\t\tUnzip ')
            utils::unzip(zipfile, exdir = dirname(genomedir))
            unlink(zipfile)
            return(invisible(genomedir))
        }
        message('\t\t\tfailed')
    }

    # Index using Bowtie
    dir.create(genomedir, showWarnings = FALSE, recursive = TRUE)
    genomefasta  <- genome_fasta(indexedgenomesdir, bsgenome)
    BSgenome::writeBSgenomeToFasta(bsgenome, genomefasta)
    Rbowtie::bowtie_build(
        genomefasta, genomedir, prefix = basename(genomedir), force = TRUE)
    return(invisible(genomedir))
    
}

#' Index targets
#' 
#' Bowtie index targets
#' @param targets   \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome  \code{\link[BSgenome]{BSgenome-class}}
#' @param outdir    string: output directory
#' @param verbose   TRUE (default) or FALSE 
#' @return invisible(targetdir)
#' @examples 
#' require(magrittr)
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' targets <- extend(bed_to_granges(bedfile, genome = 'mm10'))
#' index_targets(targets, bsgenome)
#' @export
index_targets <- function(
    targets, 
    bsgenome = getBSgenome(genome(targets)[1]),
    outdir   = OUTDIR,
    verbose  = TRUE
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
    if (verbose) cmessage(
                    '\t\t%d ranges after merging overlaps', length(targets))
    
    # Write to fasta
    targetdir   <- target_dir(outdir)
    targetfasta <- target_fasta(outdir)
    targetseqs  <- BSgenome::getSeq(bsgenome, targets)
    names(targetseqs) <- sprintf(
        '%s:%s-%s:%s',
        as.character(seqnames(targets)), start(targets), end(targets), 
        strand(targets))
    if (verbose) cmessage('\t\tWrite seqs  to %s', targetfasta)
    dir.create(dirname(targetfasta), showWarnings = FALSE, recursive = TRUE)
    Biostrings::writeXStringSet(targetseqs, targetfasta)

    # Index targets
    if (verbose) cmessage('\t\tWrite index to %s', targetdir)
    Rbowtie::bowtie_build(  targetfasta,
                            targetdir,
                            prefix = basename(targetdir),
                            force = TRUE)
    return(invisible(targetdir))
    
}


run_bowtie <- function(spacerfasta, indexdir, outfile, norc, mismatches = 2){
    
    assertive::assert_all_are_existing_files(spacerfasta)
    assertive::assert_all_are_dirs(indexdir)
    assertive::assert_is_a_bool(norc)
    assertive::assert_is_a_number(mismatches)
    assertive::assert_is_subset(mismatches, 0:3)
    
    dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
    
    Rbowtie::bowtie(
        sequences = spacerfasta,
        index     = file.path(indexdir, basename(indexdir)),
        f         = TRUE,        # fasta input
        #m         = 10000,      # ignore seqs with m+ alignments
        a         = TRUE,        # report ALL alignments
        v         = mismatches,  # up to 3 mismatches
        norc      = norc,        # no reverse complement
        outfile   = outfile,
        force     = TRUE)
}

read_bowtie_results <- function(outfile, mis){
    
    .N <- . <- readname <- NULL
    
    assertive::assert_all_are_existing_files(outfile)
    mismatch <- mismatches <- NULL
    
    dt <- data.table::fread(
            outfile,
            col.names = c(  'readname', 'strand', 'target', 'position', 
                            'readseq', 'quality', 'matches', 'mismatches'))

    pattern <- '20:[ACGT][>][ACGT]'
    dt %<>% extract(!stri_detect_regex(mismatches, pattern))

    dt[ is.na(mismatches), mismatch := 0]
    dt[!is.na(mismatches), mismatch := stringi::stri_count_fixed(
                                        mismatches, '>')]
    dt %<>% extract(mismatch<=mis)
    dt[, mismatch := factor(mismatch, c(seq(0, mis)))]

    results <-  dt %>% 
        extract( , .N, keyby = .(readname, mismatch)) %>% 
        data.table::dcast(readname ~ mismatch,  value.var='N', drop=FALSE) %>% 
        data.table::setnames(names(.)[-1], paste0('MM', names(.)[-1]))
    
    results <- cbind(results[, 1], setnafill(results[, -1], fill = 0))
    results
}


#' Match spacer sequences
#' 
#' Count matches to indexed target/genome
#' 
#' @param seqs      character vector: sequences to match against indexed ref
#' @param indexdir  string: dir containing indexed reference.
#'                  This can be an indexed genome( \code{\link{index_genome}}
#'                  It can also be indexed targets (\code{\link{index_targets}})
#' @param norc      TRUE or FALSE: whether to run bowtie also with revcompls
#'                  Generally TRUE for genome and FALSE for target matches, 
#'                  because target ranges generally include both strands.
#' @param mismatches max number of mismatches to consider
#' @param outdir    string: multicrispr output directory
#' @param verbose   TRUE (default) or FALSE
#' @return data.table 
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
#'  # indexdir <- genome_dir(indexedgenomesdir = INDEXEDGENOMESDIR, bsgenome)
#'  # match_seqs(spacers$crisprspacer, indexdir, norc=TRUE, mismatches = 1)
#'  
#' # TFBS example
#' #-------------
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' targets <- extend(bed_to_granges(bedfile, genome = 'mm10'))
#' indexdir <- index_targets(targets, bsgenome)
#' spacers <- find_spacers(targets, bsgenome)
#' seqs <- unique(paste0(spacers$crisprspacer, spacers$crisprpam))
#' match_seqs(seqs, indexdir, norc=FALSE)
#' match_seqs(seqs, indexdir, norc=FALSE, mismatches=3)
#' @export
match_seqs <- function(seqs, indexdir, norc, mismatches = 2,
    outdir = OUTDIR, verbose = TRUE
){

    # Assert
    assertive::assert_is_character(seqs)
    assertive::assert_has_no_duplicates(seqs)

    # Write reads to fasta
    reads <- Biostrings::DNAStringSet(unique(seqs))
    reads %<>% name_uniquely('read')
    readfasta <- spacer_fasta(outdir)
    dir.create(dirname(readfasta), recursive = TRUE, showWarnings = FALSE)
    if (verbose) cmessage('\t\tWrite reads to %s', readfasta)
    Biostrings::writeXStringSet(reads, readfasta)

    # Map reads and read results
    outfile <- spacer_matchfile(outdir, indexdir)
    if (verbose) cmessage('\t\tMap reads: %s', outfile)
    run_bowtie(readfasta, indexdir, outfile, norc = norc, 
                mismatches = max(1, mismatches)) # 1-mismatch offtargets 
    if (verbose) cmessage('\t\tLoad results')    # required for pam correction
    matches <- read_bowtie_results(outfile, mismatches)
    
    # Count matches
    if (verbose) cmessage('\t\tCount matches')
    readdt <- data.table(readname = names(reads), 
                        readseq   = unname(as.character(reads)))
    readdt %<>% merge(matches, by='readname', all=TRUE, sort=FALSE)
    readdt <- cbind(readdt[, c(1, 2)], 
                    data.table::setnafill(readdt[, -c(1, 2)], fill = 0))

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
#' @param outdir    string: file where to output bowtie results
#' @param pam        string (default 'NGG') pam pattern to expand
#' @param verbose    TRUE (default) or FALSE
#' @return data.table
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
#'  # indexdir <- genome_dir(indexedgenomesdir = INDEXEDGENOMESDIR, bsgenome)
#'  # match_spacers(spacers, indexdir, norc=TRUE, mismatches = 1)
#'  
#' # TFBS example
#' #-------------
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' targets <- extend(bed_to_granges(bedfile, genome = 'mm10'))
#' indexdir <- index_targets(targets, bsgenome)
#' spacers <- find_spacers(targets, bsgenome)
#' match_spacers(spacers, indexdir, norc=FALSE, mismatches = 1)
#' @export
match_spacers <- function(
    spacers, indexdir, norc, mismatches = 2,
    outdir = OUTDIR, pam = 'NGG', verbose = TRUE
){
    # Comply
    crispr <- crisprspacer <- . <- NULL
    
    # Expand pams
    if (verbose) cmessage('\t\tExpand iupac ambiguities in pam')
    spacerseqs <- unique(spacers$crisprspacer)
    pamseqs <- expand_iupac_ambiguities(pam)
    crisprdt <- data.table(
                    crisprspacer = rep(spacerseqs, each=length(pamseqs)), 
                    crisprpam    = rep(pamseqs, times=length(spacerseqs))) %>% 
                extract(, crispr := paste0(crisprspacer, pamseqs))
    crisprseqs <- unique(crisprdt$crispr)
    
    # Count matches
    matches <- match_seqs(crisprseqs, 
                            indexdir, 
                            norc       = norc, 
                            mismatches = mismatches,
                            outdir     = outdir,
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
#' @return updated spacer \code{\link[GenomicRanges]{GRanges-class}}
#' @seealso \code{\link{index_genome}}, \code{\link{index_targets}}
#' @examples
#' # TFBS example
#' #-------------
#'  bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'  bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'  targets  <- extend(bed_to_granges(bedfile, genome = 'mm10'))
#'  spacers <- find_spacers(targets, bsgenome)
#'  add_target_counts(spacers, targets, bsgenome)
#' @export
add_target_counts <- function(
    spacers, targets, bsgenome, mismatches = 2,
    pam = 'NGG', outdir = OUTDIR, verbose = TRUE
){
    # Comply
    . <- NULL
    
    # Index targets
    targetdir <- target_dir(outdir)
    dir.create(targetdir, showWarnings = FALSE, recursive = FALSE)
    index_targets(targets, bsgenome, outdir)

    # Match spacers to targets
    if (verbose) cmessage('\tAdd target counts')
    matches <- match_spacers(
                    spacers, targetdir, norc = TRUE, 
                    mismatches = mismatches, outdir = outdir,
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
#' @param spacers    spacer \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome   \code{\link[BSgenome]{BSgenome-class}}
#' @param mismatches number (default 2): max number of mismatches to consider
#' @param pam        string (default 'NGG') pam pattern to expand
#' @param outdir     dir where output is written to
#' @param indexedgenomesdir string: dir with indexed genomes
#' @param plot          FALSE (default) or TRUE
#' @param verbose       TRUE (default) or FALSE
#' @return spacer GRanges with additional mcols
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
#'  # index_genome(bsgenome)
#'  # add_genome_counts(spacers, bsgenome, mismatches=0, plot = TRUE)
#'     
#' # TFBS example
#' #-------------
#'  bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'  bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'  targets <- extend(bed_to_granges(bedfile, genome = 'mm10'))
#'  spacers <- find_spacers(targets, bsgenome)
#'  # index_genome(bsgenome)
#'  # add_genome_counts(spacers, bsgenome)
#'  # add_genome_counts(spacers, bsgenome, mismatches=3)
#' @export
add_genome_counts <- function(
    spacers, 
    bsgenome          = getBSgenome(genome(spacers)[1]),
    mismatches        = 2,
    pam               = 'NGG', 
    outdir            = OUTDIR,
    indexedgenomesdir = INDEXEDGENOMESDIR,
    plot              = FALSE,
    verbose           = TRUE
){
    . <- NULL
    
    # Add genome matches
    genomedir <- genome_dir(indexedgenomesdir, bsgenome)
    if (verbose) message('\tAdd genome counts')
    matches <- match_spacers(spacers, genomedir, norc = FALSE, 
                        mismatches = mismatches,
                        outdir = outdir, pam = pam, verbose = verbose)
    names(matches) %<>% stringi::stri_replace_first_regex('^MM', 'G')
    dt <- gr2dt(spacers) %>%
            merge(matches, by = 'crisprspacer', sort = FALSE)
    
    # Organize columns and return
    targetvars <- names(dt) %>% extract(stri_startswith_fixed(., 'target'))
    crisprvars <- c('crisprname', 'crisprspacer', 'crisprpam', 'crisprext') %>% 
                    intersect(names(dt))
    othervars <- setdiff(names(dt), c(targetvars, crisprvars))
    spacers <-  dt %>%  
                extract(, c(targetvars, crisprvars, othervars), with=FALSE) %>% 
                dt2gr(seqinfo(spacers))

    # Plot
    if (plot){
        grplot <- gr2dt(spacers) %>% 
                    extract(, if(any(.$G0==1)) .SD, by = 'targetname') %>% 
                    dt2gr(seqinfo(spacers))
        p <- plot_intervals(grplot, alpha_var = 'G0') + 
            ggplot2::scale_alpha(trans = 'reverse', range = c(0.3, 1), 
                                breaks = c(0,1,max(grplot$G0)))
        print(p)
    }
    
    # Return
    spacers
}


#' @export
#' @rdname filter_offtargets
add_specificity <- function(...){
    #.Deprecated('add_offtargets')
    add_offtargets(...)
}


#' @rdname filter_offtargets
#' @export
add_offtargets <- function(spacers, bsgenome, targets = NULL, mismatches = 2, 
    pam = 'NGG', outdir = OUTDIR, indexedgenomesdir = INDEXEDGENOMESDIR, 
    plot = TRUE, size_var = default_size_var(spacers), alpha_var = 'off', 
    verbose= TRUE
){
# First clear
    if (!has_been_indexed(bsgenome, indexedgenomesdir)) return(spacers)
    offcols <- c(paste0('G', 0:3), paste0('T', 0:3), paste0('off', 0:3), 'off')
    mcols(spacers) %<>% extract(, setdiff(names(.), offcols))
    . <- off <- NULL
# Add genome/target counts
    spacers %<>% add_genome_counts(
                    bsgenome, mismatches = mismatches, outdir = outdir,
                    indexedgenomesdir = indexedgenomesdir,
                    pam = pam, verbose = verbose)
    if (!is.null(targets)){
        spacers %<>% add_target_counts(
                        targets, bsgenome, mismatches = mismatches, 
                        outdir = outdir, pam = pam, verbose = verbose)
        for (mis in 0:mismatches){
                assert_all_are_true(mcols(spacers)[[paste0('T', mis)]] <= 
                                    mcols(spacers)[[paste0('G', mis)]])}}
# Add offtarget counts
    digits <- ceiling(log10(length(spacers)))
    if (verbose) cmessage('\tCalculate offtargets for %d spacers', 
                            length(spacers))
    spacers$off <- spacers$G0 - (if (is.null(targets)) 1 else spacers$T0)
    spacers$off0 <- spacers$off # don't switch order to keep off first
    if (!is.null(targets))  spacers$G0 <- NULL
    if (mismatches>0){
        for (mis in seq_len(mismatches)){
            Gvar <- paste0('G', mis)
            Gx <- mcols(spacers)[[Gvar]]
            Tx <- if (is.null(targets)) 0 else mcols(spacers)[[paste0('T',mis)]]
            offcounts <- Gx - Tx
            offvar <- paste0('off', mis)
            mcols(spacers)[[offvar]] <- offcounts
            spacers$off %<>% add(offcounts)
            if (is.null(targets))  mcols(spacers)[[Gvar]] <- NULL
            if (verbose) cmessage('\t       %s have no %d-mismatch offtargets', 
                format(sum(mcols(spacers)[[offvar]]==0), width = digits), mis)}}
# Plot and return
    if (plot){
        grplot <- gr2dt(spacers) %>%  # don't do that - what if no such exist?
                    #extract(, if(any(off==0)) .SD, by = 'targetname') %>% 
                    dt2gr(seqinfo(spacers))
        p <- plot_intervals(
                grplot, alpha_var = 'off', size_var = size_var)
        print(p)}
    spacers
}


#' Filter for minimal offtarget counts
#' 
#' Add mcols 'off' (offtargets), 'off0' (0-mismatch offtargets), etc.
#' 
#' For NULL targets: offtarget counts = (G0- 1) +  G1     + ...
#' Non-NULL targets: offtarget counts = (G0-T0) + (G1-T1) + ...
#' 
#' @param spacers    spacer \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome   \code{\link[BSgenome]{BSgenome-class}}
#' @param targets    NULL (default) or target 
#'                   \code{\link[GenomicRanges]{GRanges-class}}
#' @param by         string (default "targetname"): filter by this mcol
#' @param mismatches number (default 2): max number of mismatches to consider
#' @param pam        string (default 'NGG'): pam sequence
#' @param outdir     directory where output is written to
#' @param indexedgenomesdir string: dir with indexed genomes
#' @param plot       TRUE (default) or FALSE
#' @param alpha_var  string: mapped to alpha in plot
#' @param size_var   string: mapped to size in plot
#' @param verbose    TRUE (default) or FALSE
#' @param ...        to channel deprecated add_specificity to add_offtargets
#' @return  filtered spacer \code{\link[GenomicRanges]{GRanges-class}}
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
#'  # index_genome(bsgenome)
#'  spacers %>% add_offtargets(bsgenome, mismatches=0)
#'  spacers %>% filter_offtargets(bsgenome, mismatches=0)
#'  
#' # TFBS example
#' #-------------
#'  bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'  bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'  targets <- extend(bed_to_granges(bedfile, genome = 'mm10'))
#'  spacers <- find_spacers(targets, bsgenome)
#'  # index_genome(bsgenome)
#'  # spacers %>% add_offtargets(bsgenome)           # off = G - 1
#'  # spacers %>% add_offtargets(bsgenome, targets)  # off = G - T
#' @export
filter_offtargets <- function(spacers, bsgenome, targets = NULL, 
    by = 'targetname', mismatches = 2, pam = 'NGG', outdir = OUTDIR, 
    indexedgenomesdir = INDEXEDGENOMESDIR, verbose= TRUE, plot = TRUE, ...
){
    
    # Assert
    assert_is_all_of(spacers, 'GRanges')
    if (!has_been_indexed(bsgenome, indexedgenomesdir)) return(spacers)
    assert_is_subset(by, names(mcols(spacers)))

    # Add offtargets
    spacers %<>% add_offtargets(bsgenome, targets = targets, 
        mismatches = mismatches, pam = pam, outdir = outdir, 
        indexedgenomesdir = indexedgenomesdir, verbose = verbose, plot = plot, 
        ...)
    
    # Filter
    n0 <- length(spacers)
    groupby <- by
    spacers %<>% gr2dt() %>% 
                extract(, .SD[off==min(off)], by = groupby) %>% 
                dt2gr(seqinfo(spacers))
    # Message
    if (verbose) message('\tRetain ', length(spacers), '/', n0, 
                         ' ranges with minimal offtargets (per ', groupby, ')')
    
    # Return
    spacers
}

