
#=============================================================================
# bed_to_granges
#=============================================================================

#' Convert GRanges into Sequences
#' @param gr        \code{\link[GenomicRanges]{GRanges-class}}
#' @param bsgenome  \code{\link[BSgenome]{BSgenome-class}}
#' @return character vector
#' @examples 
#' bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' gr <- bed_to_granges(bedfile, bsgenome)
#' seqs(gr[1:3], bsgenome)
#' @export
seqs <- function(gr, bsgenome){
    
    # Assert
    assertive.types::assert_is_all_of(gr, 'GRanges')
    assertive.types::assert_is_all_of(bsgenome, 'BSgenome')
    
    # Do
    BSgenome::getSeq( bsgenome,
                      names  = seqnames(gr),
                      start  = start(gr),
                      end    = end(gr), 
                      strand = strand(gr), 
                      as.character = TRUE)
}
add_seqinfo <- function(gr, bsgenome){
    GenomeInfoDb::seqinfo(gr) <- GenomeInfoDb::seqinfo(bsgenome)
    gr
}


#' Get BSgenome
#' @param gr \code{\link[GenomicRanges]{GRanges-class}}
#' @return BSgenome
#' @examples 
#' bedfile  <- system.file('extdata/SRF.bed', package='multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' gr <- bed_to_granges(bedfile, bsgenome)
#' get_bsgenome(gr)
#' @export
get_bsgenome <- function(gr){
    . <- NULL
    
    assert_is_identical_to_true(is(gr, 'GRanges'))
    genome <- unique(unname(GenomeInfoDb::genome(gr)))
    assert_is_a_string(genome)
    getBSgenome(genome)
}


#' Annotate GRanges
#' @param gr \code{\link[GenomicRanges]{GRanges-class}}
#' @param txdb      \code{\link[GenomicFeatures]{TxDb-class}}
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @seealso \code{\link{genes_to_granges}} (inverse operation)
#' @examples
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' gr <- bed_to_granges(bedfile, bsgenome, plot = FALSE)
#' txdbname <- 'TxDb.Mmusculus.UCSC.mm10.knownGene'
#' txdb <- utils::getFromNamespace(txdbname, txdbname)
#' annotate_granges(gr, txdb)
#' @export
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
#' @param bedfile   file path
#' @param bsgenome \code{\link[BSgenome]{BSgenome-class}}
#' @param txdb     \code{\link[TxDb]{TXDb-class}}
#' @param addseq    logical(1)
#' @param do_order  logical(1)
#' @param plot      logical(1)
#' @param verbose   logical(1)
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @examples
#' bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' txdb <- utils::getFromNamespace('TxDb.Mmusculus.UCSC.mm10.knownGene', 
#'                                 'TxDb.Mmusculus.UCSC.mm10.knownGene')
#' gr <- bed_to_granges(bedfile, bsgenome)
#' gr <- bed_to_granges(bedfile, bsgenome, txdb)
#' @seealso \code{rtracklayer::import.bed} (documented in 
#' \code{\link[rtracklayer]{BEDFile-class}}), around which this function wraps.
#' @export
bed_to_granges <- function(
    bedfile, 
    bsgenome,
    txdb     = NULL,
    addseq   = TRUE,
    do_order = TRUE,
    plot     = TRUE, 
    verbose  = TRUE 
){
    . <- NULL
    
    # Assert
    assertive.files::assert_all_are_existing_files(bedfile)
    assertive.types::assert_is_all_of(bsgenome, 'BSgenome')
    if (!is.null(txdb))     assertive.types::assert_is_all_of(txdb, 'TxDb')
    assertive.types::assert_is_a_bool(addseq)
    assertive.types::assert_is_a_bool(do_order)
    assertive.types::assert_is_a_bool(plot)
    assertive.types::assert_is_a_bool(verbose)

    # Read
    gr <- rtracklayer::import.bed(bedfile)
    if (verbose) cmessage('\t\t%d ranges on %d chromosomes',
                    length(gr), length(unique(GenomeInfoDb::seqnames(gr))))
    
    # Set seqinfo
    assertive.sets::assert_is_subset(seqlevels(gr), seqlevels(bsgenome))
    seqlevels(gr) <- intersect(seqlevels(bsgenome), seqlevels(gr))
    seqinfo(gr) <- seqinfo(bsgenome)
    
    # Add seq
    if (addseq){
        cmessage('\t\tAdd seq')
        gr$seq <- seqs(gr, bsgenome)
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
#' @param file    Entrez Gene identifier file (one per row)
#' @param geneids Entrez Gene identifier vector
#' @param BSgenome \code{\link[BSgenome]{BSgenome-class}} or NULL. 
#'                 If specified, used to extract and add sequence to GRanges.
#' @param db      \code{\link[GenomicFeatures]{TxDb-class}} or 
#'                \code{\link[ensembldb]{EnsDb-class}}
#' @param addseq  logical(1)
#' @param plot    TRUE or FALSE
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @seealso \code{\link{annotate_granges}} (inverse operation)
#' @examples
#' 
#' # Entrez
#' #-------
#' bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#' genefile <- system.file('extdata/SRF.entrez', package='multicrispr')
#' geneids  <- as.character(read.table(genefile)[[1]])
#' db       <- utils::getFromNamespace('TxDb.Mmusculus.UCSC.mm10.knownGene',
#'                                     'TxDb.Mmusculus.UCSC.mm10.knownGene')
#' gr <- genes_to_granges(geneids, db)
#' gr <- genes_to_granges(geneids, db, bsgenome)
#' gr <- genefile_to_granges(genefile, db)
#' gr <- genefile_to_granges(genefile, db, bsgenome)
#'
#' # Ensembl
#' #--------
#' db <- EnsDb.Mmusculus.v98()
#' genefile <- system.file('extdata/SRF.ensembl', package='multicrispr')
#' geneids <- as.character(read.table(genefile)[[1]])
#' gr <- genes_to_granges(geneids, db)
#' gr <- genes_to_granges(geneids, db, bsgenome)
#' gr <- genefile_to_granges(genefile, db)
#' gr <- genefile_to_granges(genefile, db, bsgenome)
#' @export
genes_to_granges <- function(geneids, db, bsgenome = NULL, plot = TRUE){
    
    # Assert
    assertive.types::assert_is_character(geneids)
    assertive.types::assert_is_any_of(db, c('TxDb', 'EnsDb'))
    if (!is.null(bsgenome)){
        assertive.types::assert_is_any_of(bsgenome, 'BSgenome')
    }
    assertive.types::assert_is_a_bool(plot)
    
    # Convert
    gr <- GenomicFeatures::genes(db)[geneids]
    if (!is.null(bsgenome)){
        cmessage('\t\tAdd seq')
        if (seqlevelsStyle(bsgenome)[1] != seqlevelsStyle(gr)[1]){
            cmessage("\t\t\tSet seqlevelStyle(bsgenome) <- seqlevelStyle(gr)")
            seqlevelsStyle(bsgenome)[1] <- seqlevelsStyle(gr)[1]
        }
        gr$seq <- seqs(gr, bsgenome)
    }
    
    # Plot
    if (plot) plot_karyogram(gr)
    
    # Return
    gr
}


#' @rdname genes_to_granges
#' @export
genefile_to_granges <- function(file, db, bsgenome = NULL, plot = TRUE){
    assertive.files::assert_all_are_existing_files(file)
    geneids <- utils::read.table(file)[[1]] %>% as.character()
    genes_to_granges(geneids, db, bsgenome = bsgenome, plot = plot)
}


