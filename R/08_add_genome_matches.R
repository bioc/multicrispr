#' Index genome
#' 
#' Bowtie index genome
#' @param bsgenome \code{\link[BSgenome]{BSgenome-class}}
#' @param outdir string: output directory
#' @return string: directory with indexed genome
#' @examples 
#' # index_genome(BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
#' # index_genome(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
#' @export
index_genome <- function(bsgenome, outdir = '~/.multicrispr/bowtie/genome'){

    # Assert
    assert_is_all_of(bsgenome, 'BSgenome')
    
    # Create Names
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
    genomefa <- paste0(outdir, '/', bsgenome@pkgname, '.fa')
    genomedir  <- genomefa %>% substr(1, nchar(.)-3)

    # Write to fasta
    BSgenome::writeBSgenomeToFasta(bsgenome, genomefa)

    # Index genome
    subdirs <- list.dirs(genomedir, full.names = FALSE, recursive = FALSE)
    Rbowtie::bowtie_build(  genomefa,
                            genomedir,
                            prefix = basename(genomedir), 
                            force = TRUE)
    
    # Return
    return(genomedir)
}


#' Index targets
#' 
#' Bowtie index targets
#' @param targets   \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome  \code{\link[BSgenome]{BSgenome-class}}
#' @param outdir    string: output directory
#' @return
#' @examples 
#' require(magrittr)
#' bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' targets0 <- bed_to_granges(bedfile, genome = 'mm10')
#' targets  <- double_flank(targets0)
#' index_targets(targets, bsgenome, outdir = outdir)
#' @export
index_targets <- function(
    targets, bsgenome, outdir = '~/.multicrispr/bowtie/target'){

    # Write to fasta
    targetfa <- file.path(outdir, 'target.fa')
    targetdir <- targetfa %>% substr(1, nchar(.)-3)
    targets %<>% add_seq(bsgenome)
    targetseqs <- Biostrings::DNAStringSet(targets$seq)
    names(targetseqs) <- sprintf(
                          '%s:%s-%s:%s',
                          as.character(seqnames(targets)),
                          start(targets),
                          end(targets),
                          strand(targets))
    Biostrings::writeXStringSet(targetseqs, targetfa)

    # Index targets
    Rbowtie::bowtie_build(  targetfa,
                            targetdir,
                            prefix = basename(targetdir),
                            force = TRUE)
    
    # Return
    return(targetdir)
}


name_elements <- function(x, prefix = 'x'){
    names1 <- paste0(prefix, formatC(seq_along(x), 
                                digits = floor(log10(length(x))), 
                                flag = 0))
    names(x) <- names1
    x
}

run_bowtie <- function(crisprfa, indexdir, outfile, norc){
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
#' @param spacers   \code{\link[GenomicRanges]{GRanges-class}} with crispr 
#'                  sites to be matched against indexed ref
#' @param indexdir  string: dir containing indexed reference.
#'                  This can be an indexed genome( \code{\link{index_genome}}
#'                  It can also be indexed targets (\code{\link{index_targets}})
#' @param norc      TRUE or FALSE: whether to run bowtie also with revcompls
#'                  Generally TRUE for genome and FALSE for target matches, 
#'                  because target ranges generally include both strands.
#' @param outfile   string: file where to output bowtie results
#' @param outdir    string: dir where to output bowtie results
#' @param pam       string (default 'NGG') pam pattern to expand
#' @param count_vars character(3): how to name count variables
#' @param verbose   TRUE (default) or FALSE
#' @seealso \code{\link{index_genome}}, \code{\link{index_targets}}
#' @examples
#' # PE example
#' #-----------
#'     require(magrittr)
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#'     genomedir <- '~/.multicrispr/bowtie/genome/BSgenome.Hsapiens.UCSC.hg38'
#'     # index_genome(bsgenome, genomedir) # one time effort - takes few h
#'     gr <- GenomicRanges::GRanges(
#'               seqnames = c(PRNP = 'chr20:4699600',             # snp
#'                            HBB  = 'chr11:5227002',             # snp
#'                            HEXA = 'chr15:72346580-72346583',   # del
#'                            CFTR = 'chr7:117559593-117559595'), # ins
#'               strand   = c(PRNP = '+', HBB = '-', HEXA = '-', CFTR = '+'), 
#'               seqinfo  = BSgenome::seqinfo(bsgenome))
#'     spacers <- find_pe_spacers(gr, bsgenome)
#'     outdir <- '~/.multicrispr/bowtie'
#'     outfile <- file.path(outdir, 'crispr_to_ref.txt')
#'     
#'     # Match sequences
#'         match_seqs(spacers$spacer, genomedir, norc=FALSE, outfile=outfile)
#'     # Match crisprs (i.e. spacer + alternate pams)
#'         match_crisprs(spacers, genomedir, norc = FALSE, outfile = outfile)
#'     # Add "crispr to genome" matches
#'         add_genome_matches(spacers, genomedir, outdir)
#'     # Add "crispr to target" matches
#'         add_target_matches(spacers, targetdir, outdir)
#'     
#' # TFBS example
#' #-------------
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'     gdir <- '~/.multicrispr/bowtie/genome/BSgenome.Mmusculus.UCSC.mm10'
#'     # index_genome(bsgenome, gdir) # one time effort - takes few h
#'     bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'     targets0 <- bed_to_granges(bedfile, genome = 'mm10')
#'     targets  <- double_flank(targets0)
#'     tdir <- index_targets(targets, bsgenome)
#'     spacers <- find_spacers(targets, bsgenome)
#'     spacers %<>% extract(1:10) # only to keep example fast, works also else
#'     outdir <- '~/.multicrispr/bowtie'
#'     outfile <- file.path(outdir, 'crispr_to_ref.txt')
#'     
#'     # Match sequences
#'         crisprseqs <- unique(paste0(spacers$spacer, spacers$pam))
#'         match_seqs(crisprseqs, gdir, norc=FALSE, outfile=outfile)
#'     # Match crisprs (i.e. spacer + alternate pams)
#'         match_crisprs(spacers, gdir, norc = FALSE, outfile = outfile)
#'     # Add "crispr to genome" matches
#'         add_genome_matches(spacers, gdir, outdir)
#'     # Add "crispr to target" matches
#'         add_target_matches(spacers, tdir, outdir)
#' @export
match_seqs <- function(
    seqs, 
    indexdir,
    norc,
    outfile, 
    count_vars = c('MM0', 'MM1', 'MM2'),
    verbose    = TRUE
){

    # Assert
    assertive::assert_has_no_duplicates(seqs)
    
    # Write reads to fasta
    reads <- Biostrings::DNAStringSet(unique(seqs))
    reads %<>% name_elements('read')
    readfa <- file.path(outdir, 'reads.fa')
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


#' @rdname match_seqs
#' @export
match_crisprs <- function(
    spacers, indexdir, norc, outfile, pam = 'NGG', 
    count_vars = c('MM0', 'MM1', 'MM2'), verbose = TRUE){
    
    # Expand pams
    if (verbose) cmessage('\t\tExpand iupac ambiguities in pam')
    spacerseqs <- unique(spacers$spacer)
    pamseqs <- expand_iupac_ambiguities(pam)
    crisprdt <- data.table(spacer = rep(spacerseqs, each = length(pamseqs)), 
                        pam     = rep(pamseqs, times = length(spacerseqs))) %>% 
              extract(, crispr := paste0(spacer, pamseqs))
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
                               by = 'spacer')
    matches %>% setnames(c('MM0', 'MM1', 'MM2'), count_vars)
    matches[spacerseqs, on = 'spacer']
}

#' @rdname match_seqs
#' @export
add_genome_matches <- function(
    spacers, 
    indexdir,
    outdir, 
    pam        = 'NGG', 
    verbose    = TRUE
){
    outfile <- file.path(outdir, 'crispr_to_genome.txt')
    count_vars <- c('G0', 'G1', 'G2')
    matches <- match_crisprs(
                        spacers, indexdir, norc = FALSE, outfile = outfile, 
                        pam = pam, count_vars = count_vars, verbose = verbose)
    spacers %>% 
    as.data.table() %>% 
    merge(matches, by = 'spacer', sort = FALSE) %>% 
    GRanges(seqinfo = seqinfo(spacers))
}


#' @rdname match_seqs
#' @export
add_target_matches <- function(
    spacers, 
    indexdir, 
    outdir, 
    pam       = 'NGG',
    norc      = TRUE,
    verbose   = TRUE
){

    outfile <- file.path(outdir, 'crispr_to_targets.txt')
    count_vars <- c('T0', 'T1', 'T2')
    matches <- match_crisprs(
                    spacers, indexdir, norc = norc, outfile = outfile,
                    pam = pam, count_vars = count_vars, verbose = verbose)
    
    spacers %>% 
    as.data.table() %>% 
    merge(matches, by = 'spacer', sort = FALSE) %>% 
    GRanges(seqinfo = seqinfo(spacers))
}



