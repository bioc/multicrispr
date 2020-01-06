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
            col.names = c('crisprname', 'strand', 'target', 'position', 
                          'crisprseq', 'quality', 'matches', 'mismatches'))
    dt[ , mismatch := stringi::stri_count_fixed(mismatches, '>')]
    dt[ ,   list(   MM0 = sum(mismatch==0),
                    MM1 = sum(mismatch==1),
                    MM2 = sum(mismatch==2)),
            by = 'crisprname' ]
}


#' Add match counts
#' 
#' Count matches to indexed target/genome and add to GRanges
#' 
#' \code{count_matches} matches sequences against (indexed) target and genome
#' \code{add_match_counts} expands iupac amgiguities in the pam sequence, 
#' matches all resulting sequences against (indexes) target and genome, 
#' adds match counts to GRanges object, and then returns it
#' 
#' @param spacers    \code{\link[GenomicRanges]{GRanges-class}}
#' @param seqs       character vector: sequences to match to targets and genome
#' @param targetdir  string: dir with indexed target seqs (as produced by \code{\link{index_targets}})
#' @param genomedir  string: dir with indexed genome seqs (as produced by \code{\link{index_genome}}
#' @param outdir     string: dir to which to write output of mapping
#' @param verbose    TRUE (default) or FALSE
#' @examples
#' # TFBS example
#' #-------------
#' # Read and index genome
#'     require(magrittr)
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'     # genomedir <- index_genome(bsgenome) # one time operation taking few h
#'     genomedir <- '~/.multicrispr/bowtie/genome/BSgenome.Mmusculus.UCSC.mm10'
#'      
#' # Read and index targets
#'     bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#'     targets0 <- bed_to_granges(bedfile, genome = 'mm10')
#'     targets  <- double_flank(targets0)
#'     targetdir <- index_targets(targets, bsgenome)
#'     
#' # Find spacers and add match counts
#'     spacers <- find_spacers(targets, bsgenome)
#'     outdir <- '~/.multicrispr/bowtie'
#'     spacers %<>% add_match_counts(targetdir, genomedir, outdir)
#'     
#'     # This works but does not do pam expansion
#'     count_matches(unique(spacers$spacer), targetdir, genomedir, outdir)
#'     
#' # PE example
#' #-----------
#'     require(magrittr)
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#'     gr <- GenomicRanges::GRanges(
#'               seqnames = c(PRNP = 'chr20:4699600',             # snp
#'                            HBB  = 'chr11:5227002',             # snp
#'                            HEXA = 'chr15:72346580-72346583',   # del
#'                            CFTR = 'chr7:117559593-117559595'), # ins
#'               strand   = c(PRNP = '+', HBB = '-', HEXA = '-', CFTR = '+'), 
#'               seqinfo  = BSgenome::seqinfo(bsgenome))
#'     spacers <- find_pe_spacers(gr, bsgenome)
#'     spacers
#' 
#' @export
count_matches <- function(
    seqs, 
    targetdir = '~/.multicrispr/bowtie/target/target', 
    genomedir = '~/.multicrispr/bowtie/genome/genome',
    outdir    = '~/.multicrispr/bowtie', 
    verbose   = TRUE
){

    # Assert
    assertive::assert_has_no_duplicates(seqs)
    
    # Write crisprs to fasta
    crisprs <- Biostrings::DNAStringSet(unique(seqs))
    crisprs %<>% name_elements('crispr')
    crisprfa <- file.path(outdir, 'crispr.fa')
    dir.create(dirname(crisprfa), recursive = TRUE, showWarnings = FALSE)
    Biostrings::writeXStringSet(crisprs, crisprfa)

    # Align against targets
    targetout <- paste0(outdir, '/crispr_to_target.txt')
    genomeout <- paste0(outdir, '/crispr_to_genome.txt') # rc in targets!
    
    if (verbose) cmessage('\t\tMap against targets')
    run_bowtie(crisprfa, targetdir, targetout, norc = TRUE)
    if (verbose) cmessage('\t\tMap against genome')
    run_bowtie(crisprfa, genomedir, genomeout, norc = FALSE)
    
    if (verbose) cmessage('\t\tRead results')
    targetmatches <- read_bowtie_results(targetout)
    genomematches <- read_bowtie_results(genomeout)
    targetmatches %>% setnames(c('MM0', 'MM1', 'MM2'), c('T0', 'T1', 'T2'))
    genomematches %>% setnames(c('MM0', 'MM1', 'MM2'), c('G0', 'G1', 'G2'))
    
    if (verbose) cmessage('\t\tCount matches')
    crisprdt <- data.table(crisprname = names(crisprs), 
                           crisprseq = unname(as.character(crisprs)))
    crisprdt %<>% merge(targetmatches, by='crisprname', all=TRUE, sort=FALSE)
    crisprdt %<>% merge(genomematches, by='crisprname', all=TRUE, sort=FALSE)

    # Sanity check # superb - works!
    assert_all_are_true(crisprdt[, T0<=G0])
    assert_all_are_true(crisprdt[, T1<=G1])
    assert_all_are_true(crisprdt[, T2<=G2])

    # Return
    crisprdt[seqs, on = 'crisprseq']
}


#' @rdname count_matches
add_match_counts <- function(
    spacers, targetdir, genomedir, outdir, pam = 'NGG'){
    
    spacerseqs <- unique(spacers$spacer)
    pamseqs <- expand_iupac_ambiguities(pam)
    crisprdt <- data.table(spacer = rep(spacerseqs, each = length(pamseqs)), 
                        pam     = rep(pamseqs, times = length(spacerseqs))) %>% 
              extract(, crispr := paste0(spacer, pamseqs) )
    
    matchdt <- count_matches(
                    unique(crisprdt$crispr), targetdir, genomedir, outdir)
    crisprdt %<>% merge(matchdt, by.x = 'crispr', by.y = 'crisprseq')
    spacerdt <- crisprdt[, list(T0 = sum(T0), T1 = sum(T1), T2 = sum(T2), 
                                G0 = sum(G0), G1 = sum(G1), G2 = sum(G2)), 
                         by = 'spacer' ]
    spacerdt %<>% extract(spacerseqs, on = "spacer")

    spacers %<>% as.data.table() %>% 
                merge(spacerdt, by = 'spacer', sort = FALSE) %>% 
                GRanges(seqinfo = seqinfo(spacers))
    spacers
}



