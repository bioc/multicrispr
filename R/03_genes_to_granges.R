
#' Annotate GRanges
#' @param granges \code{\link[GenomicRanges]{GRanges-class}}
#' @param db      \code{\link[GenomicFeatures]{TxDb-class}} or 
#'                \code{\link[ensembldb]{EnsDb-class}}
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @seealso \code{\link{genes_to_granges}} (inverse operation)
#' @examples 
#' bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#' granges <- read_bed(bedfile, 'mm10', plot = FALSE)
#' db <- EnsDb.Mmusculus.v98()
#' annotate_granges(granges, db)
#' @export
annotate_granges <- function(granges, db){
    
    gene_id <- NULL
    
    granges %>% 
    plyranges::join_overlap_left(GenomicFeatures::genes(db)) %>% 
    data.table::as.data.table() %>% 
    extract(!is.na(gene_id) ,  
            gene_id := paste0(gene_id, collapse = ';'), 
            by = c('seqnames', 'start', 'end', 'strand')) %>% 
    unique() %>% 
    data.table::setnames('gene_id', 'is_within') %>% 
    as('GRanges')

}


#' Convert geneids into GRanges
#' @param geneids Entrezg vector
#' @param db    \code{\link[GenomicFeatures]{TxDb-class}} or 
#'                \code{\link[ensembldb]{EnsDb-class}}
#' @param plot TRUE or FALSE
#' @return \code{\link[GenomicRanges]{GRanges-class}}
#' @seealso \code{\link{annotate_granges}} (inverse operation)
#' @examples
#' # Mouse Entrez
#'     dbname <- 'TxDb.Mmusculus.UCSC.mm10.knownGene'
#'     db <- utils::getFromNamespace(dbname, dbname)
#'     geneids <- names(GenomicFeatures::genes(db))[1:100]
#'     genes_to_granges(geneids, db)
#'
#' # Human Entrez
#'     dbname <- 'TxDb.Hsapiens.UCSC.hg38.knownGene'
#'     db <- utils::getFromNamespace(dbname, dbname)
#'     geneids <- names(GenomicFeatures::genes(db))[1:100]
#'     genes_to_granges(geneids, db)
#'
#' # Mouse Ensembl
#'     db <- EnsDb.Mmusculus.v98()
#'     geneids <- names(GenomicFeatures::genes(db))[1:100]
#'     genes_to_granges(geneids, db)
#'
#' # Human Ensembl
#'     db <- EnsDb.Hsapiens.v98()
#'     geneids <- names(GenomicFeatures::genes(db))[1:100]
#'     genes_to_granges(geneids, db)
#' @export
genes_to_granges <- function(geneids, db, plot = TRUE){
    gr <- GenomicFeatures::genes(db)[geneids]
    if (plot) plot_karyogram(gr)
    gr
}


# Capitalize first letter
# @param x character vector
# @examples
# x <- 'homo'
# capitalize_first_letter(x)
# capitalize_first_letter <- function(x){
#     paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
# }
# 
# #' Parse organism
# #' @param organism 'Homo sapient'
# #' @examples
# #' parse_organism('Homo sapiens')
# #' parse_organism('homo_sapiens')
# #' parse_organism('homo.sapiens')
# parse_organism <- function(organism){
#     assertive.properties::assert_is_scalar(organism)
#     components <- stringi::stri_split_regex(organism, '[ _.]') %>% unlist()
#     assertive.properties::assert_is_of_length(components, 2)
#     
#     genus   <- components[1] %>% capitalize_first_letter()
#     species <- components[2] %>% tolower()
#     list(genus = genus, species = species)
# }
# 

# Convert ENSG vector to GRanges
# @param x ENSG vector
# @param organism 'Mus musculus', 'Homo sapiens', ...
# @param release Ensembl release. If NA (default), use current release
# @return \code{\link[GenomicRanges]{GRanges-class}}
# @export
# @examples
# x <- c('ENSMUSG00000118574', 'ENSMUSG00000118576', 'ENSMUSG00000118578')
# ensembl_to_granges(x)
# ensembl_to_granges <- function(
#     x, 
#     organism = infer_organism_from_ensembl(x), 
#     release = NA
# ){
#     
#     # Parse/Check inputs
#     assertive.types::is_character(x)
#     assertive.strings::assert_all_are_matching_fixed(x, 'ENS')
#     parsedorg <- stringi::stri_split_regex(organism, '[ . _]')
#     assertive.properties::is_of_length(parsedorg, 1)
#     orgparts <- magrittr::extract2(parsedorg, 1)
#     assertive.properties::is_of_length(orgparts, 2)
# 
#     # Make txdb object
#     txdb <- GenomicFeatures::makeTxDbFromEnsembl(
#                 organism = organism, release = release)
#     
#     # Return GRanges
#     genes(txdb)[x]
# }


# Infer organism from Ensembl gene ids
# @param x Ensembl gene Ids
# @param verbose TRUE or FALSE
# @examples 
# x <- c('ENSMUSG00000118574', 'ENSMUSG00000118576', 'ENSMUSG00000118578')
# infer_organism_from_ensembl(x)
# infer_organism_from_ensembl <- function(x, verbose = TRUE){
#     
#     # Extract prefix
#     if (verbose) message('Infer organism from ENSG identifiers ... OK')
#     assertive.types::is_character(x)
#     assertive.strings::assert_all_are_matching_fixed(x, 'ENS')
#     prefix  <-  stringi::stri_extract_first_regex(x, '(ENS[^0-9]*)G') %>% 
#                 substr(1, nchar(.)-1) %>% 
#                 unique()
#     assertive.properties::assert_is_scalar(prefix)
#     
#     # Fetch mapping table
#     url <- 'http://www.ensembl.org/info/genome/stable_ids/prefixes.html'
#     table  <-   xml2::read_html(url)      %>% 
#                 rvest::html_table()       %>% 
#                 extract2(2)               %>% 
#                 data.table::data.table()
#     
#     # Map prefix to species
#     assertive.sets::assert_is_subset(prefix, table$Prefix)
#     table[prefix, on = 'Prefix'][[2]]     %>% 
#     stringi::stri_replace_first_regex(' \\(.*\\)', '')
#     
# }


# Are all entrezg keys in db?
# @param entrezgids entrezg vector
# @param organism 'Homo sapiens', 'Mus musculus', ...
# @return TRUE or FALSE
# @examples
# entrezgids  <- c('100009600', '99889', '99982', NA_character_)
# organism <- 'Mus musculus'
# all_entrezg_are_from(entrezgids, organism)
# @export
# all_entrezg_are_from <- function(entrezgids, organism){
#     
#     x <- parse_organism(organism)
#     dbname  <- sprintf('org.%s%s.eg.db', substr(x$genus,1,1), substr(x$species,1,1))
#     db <- utils::getFromNamespace(dbname, dbname)
#     dbkeys <- AnnotationDbi::keys(db)
#     
#     entrezgids %<>% extract(!is.na(.))
#     all(entrezgids %in% dbkeys)
# }


# Infer organism from Entrezg vector
# @param x ensembl gene id vector
# @examples 
# x <- c('100009600', '99889', '99982')
# infer_organism_from_entrezg(x)
# @export
# infer_organism_from_entrezg <- function(x){
#     
#     all_organisms  <- c('Homo sapiens', 
#                         'Mus musculus', 
#                         'Danio rerio', 
#                         'Drosophila melanogaster')
#     
#     for (organism in all_organisms){
#         if (all_entrezg_are_from(x, organism)) return(organism)
#     }
#         
#     msg <- sprintf('infer_organism_from_entrezg currently limited to: %s.', 
#                     paste0(all_organisms, collapse = ', '))
#     stop(msg)
# }

