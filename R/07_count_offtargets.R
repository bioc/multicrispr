
#=============================================================================
#
#                                   PDICT_COUNT
#
#=============================================================================

#' Count matches using vcountpdict
#' 
#' @param crisprseqs character vector: sequences to match against indexed ref
#' @param referenceseqs  string: dir containing indexed referenceseqs.
#'                  This can be an indexed genome( \code{\link{index_genome}}
#'                  It can also be indexed targets (\code{\link{index_targets}})
#' @param mismatches max number of mismatches to consider
#' @param norc      TRUE or FALSE: whether to run bowtie also with revcompls
#'                  Generally TRUE for genome and FALSE for target matches, 
#'                  because target ranges generally include both strands.
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
#'  spacers <- find_primespacers(gr, bsgenome)
#'  pdict_count( spacers$crisprspacer, bsgenome,  norc=FALSE, mismatches = 0)
#'  bowtie_count(spacers$crisprspacer, index_genome(bsgenome),  norc=FALSE, 
#'              mismatches = 0)
#' @noRd
pdict_count <- function(crisprseqs, referenceseqs, mismatches, norc = FALSE,
    outdir = OUTDIR, verbose = TRUE
){
    # Assert. Comply
    assert_is_any_of(crisprseqs, c('character', 'XStringSet'))
    assert_is_any_of(referenceseqs, c('BSgenome', 'character'))
    if (is.character(referenceseqs)) assert_is_non_scalar(referenceseqs)
    assert_is_subset(mismatches, c(0,1,2,3))
    assert_is_a_bool(verbose)
    . <- count <- NULL

    # Count
    starttime <- Sys.time()
    matches <- data.table(readseq = crisprseqs)
    countfun <- if (is(referenceseqs, 'BSgenome')){     pdict_count_bsgenome
                } else if (is.character(referenceseqs)){pdict_count_character}
        
    for (i in seq(0, mismatches)){
        countvar <- paste0('MM', i)
        matches[, (countvar) := countfun(crisprseqs, referenceseqs, i)]
        if (verbose) cmessage( 
                        '\t\t\tCount %d-mismatch hits in genome      : %s',
                        i, format(signif(Sys.time() - starttime, 2)))
    }

    # Return
    matches[]
}

pdict_count_bsgenome <- function(crisprseqs, bsgenome, mismatches){
    . <- count <- NULL
    Biostrings::vcountPDict(Biostrings::DNAStringSet(crisprseqs),
                            bsgenome,
                            min.mismatch = mismatches,
                            max.mismatch = mismatches) %>% 
    data.table::as.data.table()  %>% 
    extract(, .(n = sum(count)), by ='index') %>%
    extract2('n')
}

pdict_count_character <- function(crisprseqs, targetseqs, mismatches){
    Biostrings::vcountPDict(Biostrings::DNAStringSet(crisprseqs),
                            Biostrings::DNAStringSet(targetseqs),
                            min.mismatch = mismatches,
                            max.mismatch = mismatches) %>% 
    rowSums()
}


#==============================================================================
#
#                                 BOWTIE_COUNT
#
#==============================================================================



#' Count matches using Bowtie
#' 
#' Count matches to indexed target/genome using Bowtie
#' 
#' @param crisprseqs character vector: sequences to match against indexed ref
#' @param indexdir  string: dir with Bowtie index.
#'                  This can be an indexed genome( \code{\link{index_genome}}
#'                  It can also be indexed targets (\code{\link{index_targets}})
#' @param mismatches max number of mismatches to consider
#' @param norc      TRUE or FALSE: whether to run bowtie also with revcompls
#'                  Generally TRUE for genome and FALSE for target matches, 
#'                  because target ranges generally include both strands.
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
#'  spacers <- find_primespacers(gr, bsgenome)
#'  pdict_count( spacers$crisprspacer, bsgenome,  norc=FALSE, mismatches = 0)
#'  bowtie_count(spacers$crisprspacer, index_genome(bsgenome), norc=FALSE, 
#'      mismatches = 0)
#'  
#' # TFBS example
#' #-------------
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' targets <- extend(bed_to_granges(bedfile, genome = 'mm10'))
#' indexdir <- index_targets(targets, bsgenome)
#' spacers <- find_spacers(targets[1:100], bsgenome)
#' crisprseqs <- unique(paste0(spacers$crisprspacer, spacers$crisprpam))
#' bowtie_count(crisprseqs, indexdir, norc=FALSE)
#' bowtie_count(crisprseqs, indexdir, norc=FALSE, mismatches=3)
#' @noRd
bowtie_count <- function(crisprseqs, indexdir, mismatches = 2, norc, 
    outdir = OUTDIR, verbose = TRUE
){

    # Assert
    assert_is_character(crisprseqs)
    assert_has_no_duplicates(crisprseqs)

    # Write reads to fasta
    reads <- DNAStringSet(unique(crisprseqs))
    reads %<>% name_uniquely('read')
    readfasta <- spacer_fasta(outdir)
    dir.create(dirname(readfasta), recursive = TRUE, showWarnings = FALSE)
    if (verbose) cmessage('\t\t\tWrite reads to %s', readfasta)
    writeXStringSet(reads, readfasta)

    # Map reads and read results
    outfile <- spacer_matchfile(outdir, indexdir)
    if (verbose) cmessage('\t\t\tMap reads: %s', outfile)
    run_bowtie(readfasta, indexdir, outfile, norc = norc, 
                mismatches = max(1, mismatches)) # 1-mismatch offtargets 
    if (verbose) cmessage('\t\t\tLoad results')    # required for pam correction
    matches <- read_bowtie_results(outfile, mismatches)
    
    # Count matches
    if (verbose) cmessage('\t\t\tCount matches')
    readdt <- data.table(readname = names(reads), 
                        readseq   = unname(as.character(reads)))
    readdt %<>% merge(matches, by='readname', all=TRUE, sort=FALSE)
    readdt <- cbind(readdt[, c(1, 2)], 
                    data.table::setnafill(readdt[, -c(1, 2)], fill = 0))

    # Return
    readdt[, 'readname' := NULL]
    readdt[crisprseqs, on = 'readseq']
}

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
    
    spacer_matchfile <- function(outdir = OUTDIR, bowtieindex){
        paste0(outdir, '/spacers/spacers_to_', basename(bowtieindex), '.txt')}
    
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
        res <- tryCatch(download.file(url, zipfile, method="curl"), 
                        error = function(e) 1)
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
    targets %<>% add_inverse_strand(plot = FALSE, verbose = FALSE)
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


run_bowtie <- function(spacerfasta, bowtieindex, outfile, norc, mismatches = 2){
    
    assert_all_are_existing_files(spacerfasta)
    assert_all_are_dirs(bowtieindex)
    assert_is_a_bool(norc)
    assert_is_a_number(mismatches)
    assert_is_subset(mismatches, 0:3)
    
    dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
    
    Rbowtie::bowtie(
        sequences = spacerfasta,
        index     = file.path(bowtieindex, basename(bowtieindex)),
        f         = TRUE,        # fasta input
        #m         = 10000,      # ignore seqs with m+ alignments
        a         = TRUE,        # report ALL alignments
        v         = mismatches,  # up to 3 mismatches
        norc      = norc,        # no reverse complement
        outfile   = outfile,
        force     = TRUE)
}

read_bowtie_results <- function(outfile, mis, nag_pam=TRUE, verbose=TRUE){
# Assert. Initialize.
    assert_all_are_existing_files(outfile)
    .N <- . <- readname <- mismatch <- mismatches <- NULL
# Read
    dt <- fread(outfile, col.names = c('readname', 'strand', 'target', 
                    'position', 'readseq', 'quality', 'matches', 'mismatches'))
    if (verbose)  message('\t\t\t\tRead ', nrow(dt), 
                        ' hits with max ', max(mis,1), ' mismatch(es)')
# Rm NGG "N" mismatches (avoid double counting expanded pams)
    dt %<>% extract(!stri_detect_regex(mismatches, '20:[ACGT][>][ACGT]'))
    if (verbose) message("\t\t\t\tRetain ", nrow(dt), " after removing NGG 'N'",
                        " mismatches (avoid double counting expanded pams)")
# Rm NGG "G1" mismatces (only NAG>NGG is allowed)
    dropinplus <- if (nag_pam) '21:[CT][>]G' else '21:[ACT][>]G'
    dropinmin  <- if (nag_pam) '21:[AG][>]C' else '21:[AGT][>]C'
    dt %<>% extract((strand=='+' & !stri_detect_regex(mismatches, dropinplus)) |
                    (strand=='-' & !stri_detect_regex(mismatches, dropinmin)))
        # NAG remains
        # dt[(stri_detect_regex(mismatches, '21:[ACGT]>[ACGT]'))]
    if (verbose) message("\t\t\t\tRetain ", nrow(dt), " after removing ", 
        "NGG 'G1' mismatches", 
        ifelse(nag_pam, "(except for NGG -> NAG, which is allowed)", ""))
# Rm NGG "G2" mismatches
    dropinplus <- '22:[ACT][>]G'; dropinmin  <- '22:[ATG][>]C'
    dt %<>% extract((strand=='+' & !stri_detect_regex(mismatches, dropinplus)) |
                    (strand=='-' & !stri_detect_regex(mismatches, dropinmin)))
    if (verbose) message("\t\t\t\tRetain ", nrow(dt), 
                        " after removing NGG 'G2' mismatches")
# Handle NA values
    dt[ is.na(mismatches), mismatch := 0]
# Count no of mismatches
    dt[!is.na(mismatches), mismatch := stringi::stri_count_fixed(
                                        mismatches, '>')]
# Keep only relevant number of mismatches
    dt %<>% extract(mismatch<=mis)
    if (verbose) message("\t\t\t\tRetain ", nrow(dt), 
                        " after removing more than ", mis, " mismatches")
    dt[, mismatch := factor(mismatch, c(seq(0, mis)))]
# Cast into wide format
    results <-  dt %>% 
        extract( , .N, keyby = .(readname, mismatch)) %>% 
        data.table::dcast(readname ~ mismatch,  value.var='N', drop=FALSE) %>% 
        data.table::setnames(names(.)[-1], paste0('MM', names(.)[-1]))
    results <- cbind(results[, 1], setnafill(results[, -1], fill = 0))
# Return
    results
}


explode <- function(x) unlist(strsplit(x, character(0)))
paste_dtcols <- function(x) x [, do.call(paste0, .SD) ]
expand_iupac_ambiguities <- function(x){
    paste_dtcols( 
        as.data.table( 
            lapply(Biostrings::IUPAC_CODE_MAP[explode(x)], explode)))
}


#=============================================================================
#
#                         COUNT_SPACER_MATCHES
#
#=============================================================================

#' Count spacer matches
#' 
#' Count spacer matches to referenceseqs 
#' 
#' Expands iupac amgiguities in the pam sequence.
#' Matches all resulting sequences against (indexes) target and genome.
#' Adds match counts to GRanges object, and then returns it.
#' @param spacers   spacer \code{\link[GenomicRanges]{GRanges-class}}
#' @param indexdir Bowtie index directory
#' @param mismatches number (default 2): max number of mismatches to consider
#' @param pam        string (default 'NGG') pam pattern to expand
#' @param norc      TRUE or FALSE: whether to run bowtie also with revcompls
#'                  Generally TRUE for genome and FALSE for target matches, 
#'                  because target ranges generally include both strands.
#' @param outdir    string: file where to output bowtie results
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
#'  spacers <- find_primespacers(gr, bsgenome)
#'  genomeindex <- index_genome(bsgenome)
#'  count_spacer_matches(spacers, genomeindex, norc=FALSE, mismatches=0)
#'  
#' # TFBS example
#' #-------------
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' targets <- extend(bed_to_granges(bedfile, genome = 'mm10'))
#' indexdir <- index_targets(targets, bsgenome)
#' spacers <- find_spacers(targets, bsgenome)
#' count_spacer_matches(spacers, indexdir, norc=FALSE, mismatches = 1)
#' @noRd
count_spacer_matches <- function(
    spacers, indexdir, mismatches = 2, pam = 'NGG', norc, outdir = OUTDIR, 
    verbose = TRUE
){
    # Comply
    assert_is_subset(mismatches, seq(-1,3))
    crispr <- crisprspacer <- . <- NULL
    
    # Expand pams
    if (verbose) cmessage('\t\t\tExpand iupac ambiguities in pam')
    spacerseqs <- unique(spacers$crisprspacer)
    pamseqs <- expand_iupac_ambiguities(pam)
    crisprdt <- data.table(
                    crisprspacer = rep(spacerseqs, each=length(pamseqs)), 
                    crisprpam    = rep(pamseqs, times=length(spacerseqs))) %>% 
                extract(, crispr := paste0(crisprspacer, pamseqs))
    crisprseqs <- unique(crisprdt$crispr)
    
    # Count matches
    matches <-  bowtie_count(crisprseqs, 
                        indexdir   = indexdir, 
                        mismatches = mismatches,
                        norc       = norc, 
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


#==============================================================================
#
#      count_offtargets:    add_target_matches
#                           add_genome_matches
#                           add_offtargets
#                           filter_offtargets
#
#==============================================================================


#' Add target matches
#' @param spacers     GRanges
#' @param targets     GRanges
#' @param bsgenome    BSgenome
#' @param mismatches  number
#' @param pam         string
#' @param outdir      bowtie output directory
#' @param verbose     TRUE (default) or FALSE
#' @return GRanges
#' @examples 
#' require(magrittr)
#' file <- system.file('extdata/SRF.bed', package='multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' targets0 <- bed_to_granges(file, 'mm10')
#' targets <- extend(targets0)
#' spacers <- find_spacers(targets, bsgenome, complement = FALSE, 
#'                         ontargetmethod = NULL, offtargetmethod = NULL)
#' spacers %<>% add_target_matches(targets, bsgenome)
#' @export
add_target_matches <- function(
    spacers, targets, bsgenome, mismatches = 2, pam = 'NGG', 
    outdir = OUTDIR, verbose = TRUE
){
    # Comply
    . <- NULL
    
    # Index targets
    indexdir <- target_dir(outdir)
    dir.create(indexdir, showWarnings = FALSE, recursive = FALSE)
    index_targets(targets, bsgenome, outdir)

    # Match spacers to targets
    if (verbose) cmessage('\t\tCount target cross-(mis)matches')
    matches <- count_spacer_matches(
                    spacers, indexdir, mismatches = mismatches, pam = pam, 
                    norc = TRUE, outdir = outdir, verbose = verbose)
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

#' Add genome matches
#' @param spacers            GRanges
#' @param bsgenome           BSgenome
#' @param mismatches         number
#' @param pam                string
#' @param offtargetmethod    'bowtie' or 'pdict'
#' @param outdir             bowtie output directory
#' @param indexedgenomesdir  directory with indexed genomes
#' @param verbose            TRUE (default) or FALSE
#' @return GRanges
#' @examples 
#' require(magrittr)
#' file <- system.file('extdata/SRF.bed', package='multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' targets0 <- bed_to_granges(file, 'mm10')
#' targets <- extend(targets0)
#' spacers <- find_spacers(targets, bsgenome, complement = FALSE, 
#'                         ontargetmethod = NULL, offtargetmethod = NULL)
#' spacers %<>% extract(1:100)
#' spacers %<>% add_genome_matches(bsgenome)
#' @export
add_genome_matches <- function(
    spacers, 
    bsgenome          = getBSgenome(genome(spacers)[1]),
    mismatches        = 2,
    pam               = 'NGG', 
    offtargetmethod   = c('bowtie', 'pdict')[1],
    outdir            = OUTDIR,
    indexedgenomesdir = INDEXEDGENOMESDIR,
    verbose           = TRUE
){
    assert_is_subset(offtargetmethod, c('bowtie', 'pdict'))
    if (!has_been_indexed(bsgenome, indexedgenomesdir)) return(spacers)
    if (offtargetmethod=='bowtie'){
        add_genome_matches_bowtie(
            spacers  = spacers, bsgenome = bsgenome, mismatches = mismatches,
            pam = pam, outdir = outdir, indexedgenomesdir = indexedgenomesdir, 
            verbose = verbose)
    } else {
        add_genome_matches_pdict(spacers = spacers, bsgenome = bsgenome, 
                                mismatches = mismatches, verbose = verbose)
    }
}

add_genome_matches_bowtie <- function(
    spacers, 
    bsgenome          = getBSgenome(genome(spacers)[1]),
    mismatches        = 3,
    pam               = 'NGG', 
    outdir            = OUTDIR,
    indexedgenomesdir = INDEXEDGENOMESDIR,
    verbose           = TRUE
){
    . <- NULL
    
    # Add genome matches
    indexdir <- genome_dir(indexedgenomesdir, bsgenome)
    if (verbose) message('\t\tCount genome (mis)matches')
    matches <- count_spacer_matches(
                    spacers, indexdir, mismatches = mismatches, pam = pam, 
                    norc = FALSE,  outdir = outdir, verbose = verbose)
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

    # Return
    spacers
}


add_genome_matches_pdict <- function(spacers, bsgenome, mismatches, verbose){
    . <- .N <- index <- NULL
    
    for (i in seq(0, mismatches)){
        
        # Find spacer (mis)matches
        if (verbose) message('\t\t Find ', i, '-mismatch hits for ', 
                            length(spacers), ' spacers')
        vres <- Biostrings::vmatchPDict(
            Biostrings::DNAStringSet(paste0(spacers$crisprspacer)),
            bsgenome,
            min.mismatch = i, 
            max.mismatch = i)
        seqinfo(vres) <- seqinfo(bsgenome)
        if (verbose) message('\t\t\tFound  ', length(vres), ' off', i, 
            ' hits in ', length(unique(vres$index)), ' spacer(s)')

        # Remove non-G1 matches
        vres$G1 <- BSgenome::getSeq(
            bsgenome, down_flank(vres, +2, +2), as.character = TRUE)
        vres %<>% extract(.$G1=='G' | .$G1=='A')
        if (verbose) message(
            '\t\t\tRetain ', length(vres), ' after removing non-G1 matches')
        
        # Remove non-G2 matches
        vres$G2 <- BSgenome::getSeq(
            bsgenome, down_flank(vres, +3, +3), as.character = TRUE)
        vres %<>% extract(.$G2=='G')
        if (verbose) message(
            '\t\t\tRetain ', length(vres), ' after removing non-G2 matches')
        
        # Count
        vresdt <- gr2dt(vres)[,.(MM = .N),by='index']
        Gvar <- paste0('G', i)
        data.table::setnames(vresdt, 'MM', Gvar)
        
        # Merge and return
        spacerdt <- gr2dt(spacers)
        spacerdt[, index := seq_len(.N)]
        spacerdt %<>% merge(vresdt, by = 'index', sort = FALSE, all = TRUE)
        spacerdt %>% data.table::setnafill(cols = Gvar, fill = 0)
        spacers <- dt2gr(spacerdt, seqinfo = seqinfo(spacers))
        
    }
    spacers
}



add_offtargets <- function(spacers, targets, mismatches, verbose){
    if (verbose) cmessage('\t\tCount off-targets', 
                            length(spacers))
    digits <- ceiling(log10(length(spacers)))
    spacers$off <- spacers$G0 - (if (is.null(targets)) 1 else spacers$T0)
    spacers$off0 <- spacers$off # don't switch order to keep off first
    if (verbose) cmessage('\t\t       %s are off0 free', 
            format(sum(mcols(spacers)[['off0']]==0), width = digits))
    if (is.null(targets))  spacers$G0 <- NULL
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
            if (verbose) message('\t\t       ', 
                    format(sum(mcols(spacers)[[offvar]]==0), width = digits), 
                    ' are ', offvar, ' free')}}
    return(spacers)
}

filter_offtargets <- function(
    spacers, offtargetfilterby, mismatches, verbose = TRUE
){
    if (length(offtargetfilterby)>0){
        if (verbose) message('\tFilter for best offtarget counts')
        spacerdt <- gr2dt(spacers)
        for (i in seq(0, mismatches)){
            offvar <- paste0("off", i)
            spacerdt %<>% extract(
                , .SD[get(offvar)==min(get(offvar))], by = offtargetfilterby)
            if (verbose) message('\t\tRetain ', nrow(spacerdt), ' with best ', 
                                offvar, ' counts per ', offtargetfilterby)
        }
        spacers <- dt2gr(spacerdt, seqinfo = seqinfo(spacers))
    }
    spacers
}


#' Count offtargets
#' 
#' Count (and optionally filter) offtargets
#' 
#' Add mcols 'off' (offtargets), 'off0' (0-mismatch offtargets), etc.
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
#' @param offtargetmethod     'bowtie' or 'pdict'
#' @param outdir     directory where output is written to
#' @param indexedgenomesdir string: dir with indexed genomes
#' @param plot       TRUE (default) or FALSE
#' @param verbose    TRUE (default) or FALSE
#' @param ...        passed to plot_intervals
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
#'  spacers <- find_primespacers(gr, bsgenome)
#'  # index_genome(bsgenome)
#'  count_offtargets(spacers, bsgenome, mismatches=0)
#'  
#' # TFBS example
#' #-------------
#'  bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'  bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'  targets <- extend(bed_to_granges(bedfile, genome = 'mm10'))
#'  spacers <- find_spacers(targets, bsgenome)
#'  spacers %<>% extract(1:500)
#'  # index_genome(bsgenome)
#'  # count_offtargets(spacers, bsgenome)    # off = G - 1
#'  # count_offtargets(spacers, bsgenome, targets)  # off = G - T
#' @noRd
count_offtargets <- function(spacers, bsgenome, targets = NULL, mismatches = 2, 
    pam = 'NGG', offtargetmethod = c('bowtie', 'pdict')[1], 
    offtargetfilterby = character(0), outdir = OUTDIR, 
    indexedgenomesdir = INDEXEDGENOMESDIR, verbose = TRUE, plot = TRUE, ...
){
# First clear
    if (is.null(offtargetmethod)) return(spacers)
    if (!has_been_indexed(bsgenome, indexedgenomesdir)) return(spacers)
    if (mismatches==-1) return(spacers)
    offcols <- c(paste0('G', 0:3), paste0('T', 0:3), paste0('off', 0:3), 'off')
    mcols(spacers) %<>% extract(, setdiff(names(.), offcols))
    . <- off <- NULL
# Count genome matches
    if (verbose) message('\tCount offtargets with max ', mismatches, 
                        ' mismatches using ', offtargetmethod)
    spacers %<>% add_genome_matches(
                    bsgenome, mismatches = mismatches, pam = pam, 
                    offtargetmethod = offtargetmethod, outdir = outdir,
                    indexedgenomesdir = indexedgenomesdir, verbose = verbose)
# Count target matches
    if (!is.null(targets)){
        spacers %<>% add_target_matches(
                        targets, bsgenome, mismatches = mismatches,
                        pam = pam, outdir = outdir, verbose = verbose)
        for (mis in 0:mismatches){
                assert_all_are_true(mcols(spacers)[[paste0('T', mis)]] <= 
                                    mcols(spacers)[[paste0('G', mis)]])}}
# Count and filter offtargets
    spacers %<>% add_offtargets(targets, mismatches = mismatches, 
                    verbose = verbose)
    spacers %<>% filter_offtargets(
                    offtargetfilterby = offtargetfilterby, 
                    mismatches = mismatches, verbose = verbose)
# Plot and return
    if (plot)   print(plot_intervals(spacers, ...))
    spacers
}

