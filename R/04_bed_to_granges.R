
#=============================================================================
# bed_to_granges
#=============================================================================

# Convert GRanges into Sequences
# @param gr        \code{\link[GenomicRanges]{GRanges-class}}
# @param bsgenome  \code{\link[BSgenome]{BSgenome-class}}
# @return character vector
# @examples 
# bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
# bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
# gr <- bed_to_granges(bedfile, bsgenome)
# seqs(gr[1:3], bsgenome)
# @export
# seqs <- function(gr, bsgenome){
#     
#     # Assert
#     assert_is_all_of(gr, 'GRanges')
#     assert_is_all_of(bsgenome, 'BSgenome')
#     
#     # Do
#     BSgenome::getSeq( bsgenome,
#                       names  = seqnames(gr),
#                       start  = start(gr),
#                       end    = end(gr), 
#                       strand = strand(gr), 
#                       as.character = TRUE)
# }



# Get BSgenome
# @param gr \code{\link[GenomicRanges]{GRanges-class}}
# @return BSgenome
# @examples 
# bedfile  <- system.file('extdata/SRF.bed', package='multicrispr')
# bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
# gr <- bed_to_granges(bedfile, bsgenome)
# get_bsgenome(gr)
# @export
#get_bsgenome <- function(gr){
#    . <- NULL
#    
#    assert_is_identical_to_true(is(gr, 'GRanges'))
#    genome <- unique(unname(genome(gr)))
#    assert_is_a_string(genome)
#    getBSgenome(genome)
#}


annotate_granges <- function(gr, txdb){
    
    # Assert. Import. Comply.
    assert_is_all_of(gr, 'GRanges')
    assert_is_all_of(txdb, 'TxDb')
    gene_id <- NULL

    # Align seqlevelStyle if required
    if (seqlevelsStyle(gr) != seqlevelsStyle(txdb)){
        message("Setting seqlevelsStyle(txdb) <- seqlevelsStyle(gr)")
        seqlevelsStyle(txdb) <- seqlevelsStyle(gr)
    }
    
    # Drop seqinfo (to overlap smoothly)
    txranges <- GenomicFeatures::genes(txdb)                                %>%
                as.data.table()                                             %>%
                extract(seqlevelsInUse(gr), on = 'seqnames')                %>%
                extract(, c('seqnames', 'start', 'end', 'strand', 'gene_id'),
                            with = FALSE)                                   %>%
                as('GRanges')
    
    # Overlap
    granno  <-  as.data.table(gr)                                       %>%
                as('GRanges')                                           %>%
                plyranges::join_overlap_left(txranges)                  %>%
                as.data.table()                                         %>%
                extract(!is.na(gene_id) ,  
                        gene_id := paste0(gene_id, collapse = ';'), 
                        by = c('seqnames', 'start', 'end', 'strand'))   %>%
                unique()                                                %>%
                as('GRanges')
    seqlevels(granno) %<>% setdiff('.') # patch plyranges bug
    seqinfo(granno) <- seqinfo(gr)
    granno
}



#' Read bedfile into GRanges
#' @param bedfile   file path
#' @param genome    string: UCSC genome name (e.g. 'mm10')
#' @param txdb      NULL (default) or \code{\link[GenomicFeatures]{TxDb-class}}
#'                  (used for gene annotation)
#' @param do_order  TRUE (default) or FALSE: order on seqnames and star?
#' @param plot      TRUE (default) or FALSE: plot karyogram?
#' @param verbose   TRUE (default) or FALSE
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @seealso \code{\link{char_to_granges}}, \code{\link{genes_to_granges}}
#' @examples
#' bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' (gr <- bed_to_granges(bedfile, genome='mm10'))
#' @export
bed_to_granges <- function(
    bedfile,
    genome,
    txdb       = NULL,
    do_order   = TRUE,
    plot       = TRUE, 
    verbose    = TRUE 
){
    . <- NULL
    
    # Assert
    assert_all_are_existing_files(bedfile)
    if (!is.null(txdb)) assert_is_all_of(txdb, 'TxDb')
    assert_is_a_bool(do_order)
    assert_is_a_bool(plot)
    assert_is_a_bool(verbose)

    # Read
    if (verbose) cmessage('\tRead %s into GRanges', basename(bedfile))
    gr <- rtracklayer::import.bed(bedfile, genome = genome)
    if (verbose) cmessage('\t\t%d ranges on %d chromosomes',
                    length(gr), length(unique(seqnames(gr))))
    
    # Annotate
    if (!is.null(txdb)){
        cmessage('\t\tAnnotate with txdb')
        gr %<>% annotate_granges(txdb)
    }
    
    # Plot
    genome1 <- unique(genome(gr))
    assertive.properties::assert_is_scalar(genome1)
    title <- paste0(genome1, ': ', basename(bedfile))
    if (plot) plot_karyogram(gr, title)
    
    # Order
    if (do_order)  gr %<>% sort(ignore.strand = TRUE)
                    #%<>% extract( order(seqnames(.), start(.)))

    # Record    
    names(gr) <- gr$targetname  <- make_unique_names(gr, 'T')
    gr$targetstart <- GenomicRanges::start(gr)
    gr$targetend   <- GenomicRanges::end(gr)

    # Return
    gr
}


#' Convert character vector into GRanges
#' @param x character vector
#' @param bsgenome \code{\link[BSgenome]{BSgenome-class}}
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @examples 
#' require(magrittr)
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#' x <- c(PRNP  = 'chr20:4699600:+',            # snp
#'        HBB  = 'chr11:5227002:-',            # snp
#'        HEXA = 'chr15:72346580-72346583:-',  # del
#'        CFTR = 'chr7:117559593-117559595:+') # ins
#' gr <- char_to_granges(x, bsgenome)
#' plot_intervals(gr, facet_var = c('targetname', 'seqnames'))
#' @seealso \code{\link{bed_to_granges}}, \code{\link{genes_to_granges}}
#' @export
char_to_granges <- function(x, bsgenome){
    gr <- GenomicRanges::GRanges(x, seqinfo  = BSgenome::seqinfo(bsgenome))
    names(gr) <- gr$targetname  <- make_unique_names(gr, 'T')
    gr$targetstart <- GenomicRanges::start(gr)
    gr$targetend   <- GenomicRanges::end(gr)
    gr
}



#' Convert geneids into GRanges
#' @param file       Gene identifier file (one per row)
#' @param geneids    Gene identifier vector
#' @param complement TRUE (default) or FALSE: add complementary strand?
#' @param txdb       \code{\link[GenomicFeatures]{TxDb-class}} or 
#'                   \code{\link[ensembldb]{EnsDb-class}}
#' @param plot       TRUE (default) or FALSE
#' @param verbose    TRUE (default) or FALSE
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @seealso \code{\link{char_to_granges}}, \code{\link{bed_to_granges}}
#' @examples
#' # Entrez
#' #-------
#'     genefile <- system.file('extdata/SRF.entrez', package='multicrispr')
#'     geneids  <- as.character(read.table(genefile)[[1]])
#'     txdb     <- getFromNamespace('TxDb.Mmusculus.UCSC.mm10.knownGene',
#'                              'TxDb.Mmusculus.UCSC.mm10.knownGene')
#'     (gr <- genes_to_granges(geneids, txdb))
#'     (gr <- genefile_to_granges(genefile, txdb))
#'
#' # Ensembl
#' #--------
#'     txdb <- EnsDb.Mmusculus.v98()
#'     genefile <- system.file('extdata/SRF.ensembl', package='multicrispr')
#'     geneids <- as.character(read.table(genefile)[[1]])
#'     (gr <- genes_to_granges(geneids, txdb))
#'     (gr <- genefile_to_granges(genefile, txdb))
#' @export
genes_to_granges <- function(
    geneids, 
    txdb, 
    complement = TRUE, 
    plot       = TRUE, 
    verbose    = TRUE
){
    
    # Assert
    assert_is_character(geneids)
    assert_is_any_of(txdb, c('TxDb', 'EnsDb'))
    assert_is_a_bool(complement)
    assert_is_a_bool(plot)
    
    # Convert
    gr <- GenomicFeatures::genes(txdb)[geneids]
    if (verbose)  cmessage('\t\tConvert %d genes to %d GRanges', 
                            length(geneids), length(gr))
    
    # Add complementary strand
    if (complement){
        gr %<>% add_inverse_strand(plot = FALSE, verbose = verbose)
    }

    # Plot
    if (plot) plot_karyogram(gr)
    
    # Record    
    names(gr) <- gr$targetname  <- make_unique_names(gr, 'T')
    gr$targetstart <- GenomicRanges::start(gr)
    gr$targetend   <- GenomicRanges::end(gr)

    # Return
    gr
}


#' @rdname genes_to_granges
#' @export
genefile_to_granges <- function(file, txdb, complement = TRUE, plot = TRUE){
    assert_all_are_existing_files(file)
    geneids <- read.table(file)[[1]] %>% as.character()
    genes_to_granges(geneids, txdb, complement = complement, plot = plot)
}


