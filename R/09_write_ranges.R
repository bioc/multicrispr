#' Write GRanges to file
#' @param gr GRanges
#' @param file 
#' @param verbose TRUE (default) or FALSE
#' @return NULL
#' @examples 
#' # Find PE spacers for 4 clinically relevant loci (Anzalone et al, 2019)
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#'     gr <- char_to_granges(c(
#'         PRNP = 'chr20:4699600:+',             # snp: prion disease
#'         HBB  = 'chr11:5227002:-',             # snp: sickle cell anemia
#'         HEXA = 'chr15:72346580-72346583:-',   # del: tay sachs disease
#'         CFTR = 'chr7:117559593-117559595:+'), # ins: cystic fibrosis
#'         bsgenome)
#'     spacers <- find_pe_spacers(gr, bsgenome)
#'     targetfile <- file.path(tempdir(), 'targets.txt')
#'     spacerfile <- file.path(tempdir(), 'spacers.txt')
#'     write_ranges(gr,      targetfile)
#'     write_ranges(spacers, spacerfile)
#'     read_ranges(targetfile)
#' @export
write_ranges <- function(gr, file, verbose=TRUE){
    assert_all_are_dirs(dirname(file))
    assert_is_all_of(gr, 'GRanges')
    
    fwrite(gr2dt(gr), file, sep = '\t')
    if (verbose) message('Write to ', file)
}

#' @rdname write_ranges
#' @export
read_ranges <- function(file){
    data.table::fread(file)
}