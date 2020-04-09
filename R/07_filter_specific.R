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

crisprseq_matchfile <- function(outdir = OUTDIR, indexdir){
        paste0(outdir, '/crisprseqs/crisprseqs_to_', basename(indexdir), '.txt')}

crispr_fasta <- function(outdir = OUTDIR){
        paste0(outdir, '/crisprseqs.fa') }

#' Has been indexed?
#' @param bsgenome BSgenome
#' @param indexedgenomesdir directory with indexed genomes
#' @examples 
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#' has_been_indexed(bsgenome)
#' @export
has_been_indexed <- function(bsgenome, indexedgenomesdir = INDEXEDGENOMESDIR){
    dir.exists(genome_dir(bsgenome = bsgenome))
}


#' Index genome
#' 
#' Bowtie index genome
#' @param bsgenome \code{\link[BSgenome]{BSgenome-class}}
#' @param indexedgenomesdir  string: directory with bowtie-indexed genome
#' @param overwrite logical(1)
#' @return invisible(genomdir)
#' @examples
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' #index_genome(bsgenome)
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
#' #index_genome(bsgenome)
#' @export
index_genome <- function(
    bsgenome, indexedgenomesdir=INDEXEDGENOMESDIR, overwrite=FALSE
){
    # Assert
    assert_is_all_of(bsgenome, 'BSgenome')
    
    # Return if subdir exists
    genomedir    <- genome_dir(  indexedgenomesdir, bsgenome)
    genomefasta  <- genome_fasta(indexedgenomesdir, bsgenome)
    if (!overwrite & 
        dir.exists(genomedir) & length(list.files(genomedir))!=0){
        cmessage('%s already contains index - set overwrite=TRUE to overwrite', 
                genomedir)
        return(invisible(indexedgenomesdir))
    }
    
    # Create subdir and write fastafile and genomeindex
    dir.create(genomedir, showWarnings = FALSE, recursive = TRUE)
    BSgenome::writeBSgenomeToFasta(bsgenome, genomefasta)
    Rbowtie::bowtie_build(
        genomefasta, genomedir, prefix = basename(genomedir), force = TRUE)
    return(invisible(indexedgenomesdir))
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

rm_pam_mismatches <- function(dt, pam){
    mismatches <- NULL
    if (pam=='NGG'){
        pattern <- '20:[ACGT][>][ACGT]'
        dt %<>% extract(!stri_detect_regex(mismatches, pattern))
    } else {
        message('\t\tOfftarget analysis requires additional work ', 
                'for pam patterns other than NGG ')
    }
    return(dt)
}

run_bowtie <- function(
    crisprfasta, indexdir, outfile, norc, mismatches = 2, pam = 'NGG'
){
    # Assert. Comply.
    assert_all_are_existing_files(c(crisprfasta))
    assert_all_are_dirs(indexdir)
    assert_is_a_bool(norc)
    assert_is_a_number(mismatches)
    assert_is_subset(mismatches, 1:3) # 1 for alternate pams
    assert_is_a_string(pam)
    .N <- . <- spacername <- crisprname <- mismatch <- NULL
    
    # Run bowtie
    dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
    bowtie( sequences = crisprfasta,
            index     = file.path(indexdir, basename(indexdir)),
            f         = TRUE,        # fasta input
            #m         = 10000,      # ignore seqs with m+ alignments
            a         = TRUE,        # report ALL alignments
            v         = mismatches,  # up to 3 mismatches
            norc      = norc,        # no reverse complement
            outfile   = outfile,
            force     = TRUE)

    # Read results
    dt <- fread(outfile,
                col.names = c('crisprname', 'strand', 'target', 'position', 
                            'crisprseq', 'quality', 'matches', 'mismatches'))
    dt[, spacername := crisprname %>% tstrsplit('.', fixed=TRUE) %>% extract(1)]
    dt[, pam        := crisprname %>% tstrsplit('.', fixed=TRUE) %>% extract(2)]
    
    # Count matches
    dt %<>% rm_pam_mismatches(pam) # dont double count pams!
    dt[ is.na(mismatches), mismatch := 0]
    dt[!is.na(mismatches), mismatch := stri_count_fixed(mismatches, '>')]
    results <-  dt[ , .N, keyby = .(spacername, mismatch)] %>% 
                data.table::dcast(spacername ~ mismatch,  value.var = 'N') %>% 
                setnames(names(.)[-1], paste0('MM', names(.)[-1]))
    results <- cbind(results[, 1], setnafill(results[, -1], fill = 0))
    
    # Return
    results
}


explode <- function(x) unlist(strsplit(x, character(0)))
paste_dtcols <- function(x) x [, do.call(paste0, .SD) ]
expand_pam <- function(pam){
    assert_is_a_string(pam)
    paste_dtcols( 
        as.data.table( 
            lapply(Biostrings::IUPAC_CODE_MAP[explode(pam)], explode)))
}

factorcode <- function(x, prefix = 's'){
    . <- NULL
    x %>% 
    factor(unique(.)) %>% 
    as.integer() %>% 
    formatC(digits = floor(log10(max(.))), flag=0) %>% paste0(prefix, .)
}

#' @rdname add_genome_counts
#' @export
add_match_counts <- function(spacers, indexdir, norc, mismatches = 2,
    outdir = OUTDIR, pam = 'NGG', verbose = TRUE
){
    # Assert. Comply
    assert_is_all_of(spacers, 'GRanges')
    assert_all_are_dirs(c(indexdir, outdir))
    assert_is_a_bool(verbose)
    crispr <- crisprspacer <- crisprpam <- crisprname <- spacername <- . <- NULL

    # Get spacers. Expands pams. Construct crisprseqs.    
    spacerdt <- data.table(crisprspacer = unique(spacers$crisprspacer))
    spacerdt %>% extract(, spacername := factorcode(crisprspacer))
    if (verbose) cmessage('\t\tExpand iupac ambiguities in pam')
    crisprdt <- data.table(expand.grid(
                                crisprpam = expand_pam(pam), 
                                crisprspacer = spacerdt$crisprspacer))[,2:1]
    crisprdt %<>% merge(spacerdt, by = 'crisprspacer', sort = FALSE)
    crisprdt %>% extract(, crispr     := paste0(crisprspacer,    crisprpam))
    crisprdt %>% extract(, crisprname := paste0(spacername, '.', crisprpam))
    
    # Write crisprseqs to fasta
    crisprseqs <- Biostrings::DNAStringSet(crisprdt$crispr)
    names(crisprseqs) <- crisprdt$crisprname
    crisprfasta <- crispr_fasta(outdir)
    dir.create(dirname(crisprfasta), recursive = TRUE, showWarnings = FALSE)
    if (verbose) cmessage('\t\tWrite crisprseqs to %s', crisprfasta)
    Biostrings::writeXStringSet(crisprseqs, crisprfasta)

    # Run bowtie and count matches
    outfile <- crisprseq_matchfile(outdir, indexdir)
    if (verbose) cmessage('\t\tBowtie crisprseqs to %s: %s', 
                        basename(indexdir), outfile)
    matches <- run_bowtie(
        crisprfasta, indexdir, outfile, norc = norc, mismatches = mismatches, pam = pam)
    spacerdt %<>% merge(matches, by='spacername', all=TRUE, sort=FALSE)
    mcols(spacers) %<>% merge(spacerdt[, -1], by = 'crisprspacer', all = TRUE,
                            sort = FALSE)

    # Organize columns and return
    targetvars <- names(mcols(spacers))[stri_startswith_fixed(., 'target')]
    crisprvars <- c('crisprname', 'crisprspacer', 'crisprpam', 'crisprext') %>% 
                    intersect(names(mcols(spacers)))
    othervars <- setdiff(names(mcols(spacers)), c(targetvars, crisprvars))
    mcols(spacers) %<>% extract(, c(targetvars, crisprvars, othervars))
    spacers
}


#' @rdname add_genome_counts
#' @export
add_target_counts <- function(
    spacers, targets, bsgenome, mismatches = 2,
    pam = 'NGG', outdir = OUTDIR, verbose = TRUE
){
    
    # Index targets
    targetdir <- index_targets(targets, bsgenome, outdir = outdir)

    # Match spacers to targets
    if (verbose) cmessage('\tAdd target match counts')
    spacers %<>% add_match_counts(
                    targetdir, norc = TRUE, 
                    mismatches = mismatches, outdir = outdir,
                    pam = pam, verbose = verbose)
    names(mcols(spacers)) %<>% stringi::stri_replace_first_regex('^MM', 'T')
    spacers
}

#' Add genome/targetset match counts
#' 
#' Count spacer matches to targetset/genome and add to GRanges
#' 
#' Expands iupac amgiguities in the pam sequence.
#' Matches all resulting sequences against (indexes) target and genome.
#' Adds match counts to GRanges object, and then returns it.
#' @param spacers    spacer \code{\link[GenomicRanges]{GRanges-class}}
#' @param targets  target \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome   \code{\link[BSgenome]{BSgenome-class}}
#' @param mismatches number (default 2): max number of mismatches to consider
#' @param pam        string (default 'NGG') pam pattern to expand
#' @param outdir     dir where output is written to
#' @param indexdir   string: dir with indexed reference
#' @param indexedgenomesdir string: dir with indexed genomes
#' @param norc       TRUE or FALSE: whether to search reverse complement of 
#'                   reference seqs
#' @param verbose       TRUE (default) or FALSE
#' @return spacer GRanges with additional mcols
#' @seealso \code{\link{index_genome}}, \code{\link{index_targets}}
#' @examples
#' # TFBS example
#' #-------------
#'  bs <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'  bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'  targets  <- extend(bed_to_granges(bedfile, genome = 'mm10'))
#'  spacers <- find_spacers(targets, bs)
#'  add_target_counts(spacers, targets, bs)
#'  add_match_counts( spacers, index_targets(targets, bs), norc=FALSE)
#'  # add_genome_counts(spacers, bs, indexedgenomesdir=index_genome(bs))
#'
#' # PE example
#' #-----------
#'  require(magrittr)
#'  bs <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#'  gr <- char_to_granges(c(PRNP = 'chr20:4699600:+',             # snp
#'                          HBB  = 'chr11:5227002:-',             # snp
#'                          HEXA = 'chr15:72346580-72346583:-',   # del
#'                          CFTR = 'chr7:117559593-117559595:+'), # ins
#'                        bs)
#'  spacers <- find_pe_spacers(gr, bs)
#'  # add_genome_counts(spacers, bs, indexedgenomesdir = index_genome(bs))
#' @export
add_genome_counts <- function(
    spacers, 
    bsgenome          = getBSgenome(genome(spacers)[1]),
    mismatches        = 2,
    pam               = 'NGG', 
    outdir            = OUTDIR,
    indexedgenomesdir = INDEXEDGENOMESDIR,
    verbose           = TRUE
){
    genomedir <- genome_dir(indexedgenomesdir, bsgenome)
    if (verbose) message('\tAdd genome match counts')
    spacers %<>% add_match_counts(genomedir, norc = FALSE, 
                        mismatches = mismatches, outdir = outdir, 
                        pam = pam, verbose = verbose)
    names(mcols(spacers)) %<>% stringi::stri_replace_first_regex('^MM', 'G')
    spacers
}


#' @rdname filter_target_specific
#' @export
add_specificity <- function(
    spacers, targets, bsgenome, mismatches = 2, pam = 'NGG', outdir = OUTDIR,
    indexedgenomesdir = INDEXEDGENOMESDIR, plot = TRUE, verbose= TRUE
){
    # Comply
    specific <- NULL
    
    # Add target/genome matches
    spacers %<>% add_target_counts(
                    targets, bsgenome, mismatches = mismatches, outdir = outdir,
                    pam = pam, verbose = verbose)
    spacers %<>% add_genome_counts(
                    bsgenome, mismatches = mismatches, outdir = outdir,
                    indexedgenomesdir = indexedgenomesdir,
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
        idx %<>% and(   mcols(spacers)[[paste0('T', mis)]] == 
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
#' @param spacers    spacer \code{\link[GenomicRanges]{GRanges-class}}
#' @param targets    target \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome   \code{\link[BSgenome]{BSgenome-class}}
#' @param mismatches number (default 2): max number of mismatches to consider
#' @param pam        string (default 'NGG'): pam sequence
#' @param outdir     directory where output is written to
#' @param indexedgenomesdir string: dir with indexed genomes
#' @param plot       TRUE (default) or FALSE
#' @param verbose    TRUE (default) or FALSE
#' @return  filtered spacer \code{\link[GenomicRanges]{GRanges-class}}
#' @examples
#' # TFBS example
#' #-------------
#'  require(magrittr)
#'  bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'  bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'  targets <- extend(bed_to_granges(bedfile, genome = 'mm10'))
#'  spacers <- find_spacers(targets, bsgenome)
#'  # index_genome(bsgenome)
#'  # spacers %<>% add_specificity(targets, bsgenome)
#'  # spacers %>% filter_target_specific(targets, bsgenome)
#' @export
filter_target_specific <- function(
    spacers,
    targets,
    bsgenome          = getBSgenome(genome(spacers)[1]),
    mismatches        = 2,
    pam               = 'NGG',
    outdir            = OUTDIR, 
    indexedgenomesdir = INDEXEDGENOMESDIR,
    plot              = TRUE,
    verbose           = TRUE
){
    # Add specificty info
    spacers %<>% add_specificity(
                    targets       = targets, 
                    bsgenome      = bsgenome, 
                    mismatches    = mismatches,
                    pam           = pam, 
                    outdir        = outdir, 
                    indexedgenomesdir = indexedgenomesdir,
                    plot          = plot,
                    verbose       = verbose)
    
    # Subset
    specific <- NULL
    spacers %<>% subset(specific==TRUE)
    spacers$specific <- FALSE
    
    # Return
    return(spacers)
}



#' Filter for specific prime editing spacers
#' 
#' Filters spacers which are specific for prime editing site
#' 
#' @param spacers           spacer \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome          \code{\link[BSgenome]{BSgenome-class}}
#' @param outdir            directory where output is written to
#' @param pam               string (default 'NGG'): pam sequence
#' @param outdir            string: dir to which output is written
#' @param indexedgenomesdir dir with bowtie-indexed genomes
#' @param verbose           TRUE (default) or FALSE
#' @return filtered spacer \code{\link[GenomicRanges]{GRanges-class}}
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
#'  # index_genome(bsgenome) # one time effort - takes few h
#'  # filter_prime_specific(spacers, bsgenome)
#' @export
filter_prime_specific <- function(
    spacers, 
    bsgenome          = getBSgenome(genome(spacers)[1]),
    pam               = 'NGG', 
    outdir            = OUTDIR, 
    indexedgenomesdir = INDEXEDGENOMESDIR,
    verbose           = TRUE
){
    # Add genome matches
    spacers %<>% add_genome_counts(
                    bsgenome, mismatches = 1, outdir = outdir, 
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

