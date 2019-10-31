
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
#     assertive.types::assert_is_all_of(gr, 'GRanges')
#     assertive.types::assert_is_all_of(bsgenome, 'BSgenome')
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
#    genome <- unique(unname(GenomeInfoDb::genome(gr)))
#    assert_is_a_string(genome)
#    getBSgenome(genome)
#}


add_inverse_strand <- function(gr, plot = TRUE, verbose = TRUE){
    complements <- invertStrand(gr)
    newranges <- c(gr, complements)
    txt <- sprintf('\t\t%d ranges after adding inverse strands',
                    length(newranges))
    if (plot){
        plot_intervals(
            GRangesList(original = gr, complements = complements),
            title = txt)
    }
    if (verbose) cmessage(txt)
    newranges
}


annotate_granges <- function(gr, txdb){
    
    # Assert
    assertive.types::assert_is_all_of(gr, 'GRanges')
    assertive.types::assert_is_all_of(txdb, 'TxDb')
    gene_id <- NULL

    # Align seqlevelStyle if required    
    if (seqlevelsStyle(gr) != seqlevelsStyle(txdb)){
        message("Setting seqlevelsStyle(txdb) <- seqlevelsStyle(gr)")
        seqlevelsStyle(txdb) <- seqlevelsStyle(gr)
    }
    
    # Drop seqinfo (to overlap smoothly)
    txranges <- GenomicFeatures::genes(txdb)                                %>%
                data.table::as.data.table()                                 %>%
                extract(seqlevelsInUse(gr), on = 'seqnames')                %>%
                extract(, c('seqnames', 'start', 'end', 'strand', 'gene_id'),
                            with = FALSE)                                   %>%
                as('GRanges')
    
    # Overlap
    granno <- data.table::as.data.table(gr)                         %>%
            as('GRanges')                                           %>%
            plyranges::join_overlap_left(txranges)                  %>%
            data.table::as.data.table()                             %>%
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
#' 
#' @param bedfile    file path
#' @param genome     string: UCSC genome name (e.g. 'mm10')
#' @param complement TRUE (default) or FALSE: add complementary strand too?
#' @param txdb       NULL (default) or \code{\link[GenomicFeatures]{TxDb-class}} (for gene annotation)
#' @param do_order   TRUE (default) or FALSE: order on seqnames and star?
#' @param plot       TRUE (default) or FALSE: plot karyogram?
#' @param verbose    TRUE (default) or FALSE
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @examples
#' bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' gr <- bed_to_granges(bedfile, 'mm10')
#' @seealso \code{rtracklayer::import.bed} (documented in 
#' \code{\link[rtracklayer]{BEDFile-class}}), around which this function wraps.
#' @export
bed_to_granges <- function(
    bedfile,
    genome,
    complement = TRUE,
    txdb       = NULL,
    do_order   = TRUE,
    plot       = TRUE, 
    verbose    = TRUE 
){
    . <- NULL
    
    # Assert
    assertive.files::assert_all_are_existing_files(bedfile)
    assertive.types::assert_is_a_bool(complement)
    if (!is.null(txdb)) assertive.types::assert_is_all_of(txdb, 'TxDb')
    assertive.types::assert_is_a_bool(do_order)
    assertive.types::assert_is_a_bool(plot)
    assertive.types::assert_is_a_bool(verbose)

    # Read
    if (verbose) cmessage('\tRead %s into GRanges', basename(bedfile))
    gr <- rtracklayer::import.bed(bedfile, genome = genome)
    if (verbose) cmessage('\t\t%d ranges on %d chromosomes',
                    length(gr), length(unique(GenomeInfoDb::seqnames(gr))))
    
    # Add complementary strand
    if (complement){
        gr %<>% add_inverse_strand(plot = FALSE, verbose = verbose)
    }
    
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
    if (do_order)
        gr %<>% extract(order(seqnames(.), start(.)))
    
    # Return
    gr
}


#=============================================================================
# genes_to_granges
#=============================================================================

#' Convert geneids into GRanges
#' @param file       Gene identifier file (one per row)
#' @param geneids    Gene identifier vector
#' @param complement TRUE (default) or FALSE: add complementary strand?
#' @param txdb       \code{\link[GenomicFeatures]{TxDb-class}} or 
#'                   \code{\link[ensembldb]{EnsDb-class}}
#' @param plot       TRUE (default) or FALSE
#' @param verbose    TRUE (default) or FALSE
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @examples
#' # Entrez
#' #-------
#' genefile <- system.file('extdata/SRF.entrez', package='multicrispr')
#' geneids  <- as.character(read.table(genefile)[[1]])
#' txdb     <- utils::getFromNamespace('TxDb.Mmusculus.UCSC.mm10.knownGene',
#'                                     'TxDb.Mmusculus.UCSC.mm10.knownGene')
#' gr <- genes_to_granges(geneids, txdb)
#' gr <- genefile_to_granges(genefile, txdb)
#'
#' # Ensembl
#' #--------
#' txdb <- EnsDb.Mmusculus.v98()
#' genefile <- system.file('extdata/SRF.ensembl', package='multicrispr')
#' geneids <- as.character(read.table(genefile)[[1]])
#' gr <- genes_to_granges(geneids, txdb)
#' gr <- genefile_to_granges(genefile, txdb)
#' @export
genes_to_granges <- function(
    geneids, 
    txdb, 
    complement = TRUE, 
    plot       = TRUE, 
    verbose    = TRUE
){
    
    # Assert
    assertive.types::assert_is_character(geneids)
    assertive.types::assert_is_any_of(txdb, c('TxDb', 'EnsDb'))
    assertive.types::assert_is_a_bool(complement)
    assertive.types::assert_is_a_bool(plot)
    
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
    
    # Return
    gr
}


#' @rdname genes_to_granges
#' @export
genefile_to_granges <- function(file, txdb, complement = TRUE, plot = TRUE){
    assertive.files::assert_all_are_existing_files(file)
    geneids <- utils::read.table(file)[[1]] %>% as.character()
    genes_to_granges(geneids, txdb, complement = complement, plot = plot)
}


#' Add sequence to GRanges
#' @param gr        \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome  \code{\link[BSgenome]{BSgenome-class}}
#' @param verbose   logical(1)
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @examples 
#' bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' gr <- bed_to_granges(bedfile, 'mm10')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' gr <- add_seq(gr, bsgenome)
#' gr
#' @export
add_seq <- function(gr, bsgenome, verbose = TRUE){
    
    # Assert
    assertive.types::assert_is_all_of(gr, 'GRanges')
    assertive.types::assert_is_all_of(bsgenome, 'BSgenome')
    
    # Message
    if (verbose)  cmessage('\tAdd seq')
    
    # Align seqlevelsStyle if required
    if (seqlevelsStyle(bsgenome)[1] != seqlevelsStyle(gr)[1]){
        cmessage("\t\t\tSet seqlevelStyle(bsgenome) <- seqlevelStyle(gr)")
        seqlevelsStyle(bsgenome)[1] <- seqlevelsStyle(gr)[1]
    }
    
    # Add seq
    gr$seq <- unname(BSgenome::getSeq(
                bsgenome,
                names        = seqnames(gr),
                start        = start(gr),
                end          = end(gr), 
                strand       = strand(gr), 
                as.character = TRUE))
    
    # Return
    gr
}

