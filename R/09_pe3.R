

#' Find nicking spacers
#' @examples
#' bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#' gr <- char_to_granges(c(
#'         PRNP = 'chr20:4699600:+',             # snp: prion disease
#'         HBB  = 'chr11:5227002:-',             # snp: sickle cell anemia
#'         HEXA = 'chr15:72346580-72346583:-',   # del: tay sachs disease
#'         CFTR = 'chr7:117559593-117559595:+'), # ins: cystic fibrosis
#'         bsgenome)
#' pespacers <- find_pe_spacers(gr, bsgenome)
#' # index_genome(bsgenome)
#' # pespacers %<>% add_offtargets(bsgenome, mismatches=0)
#' find_spacers(extend_for_pe(gr), bsgenome, complement=FALSE)
find_nicking_spacers <- function(pespacers, bsgenome, plot = TRUE){
    
    pespacers$extension     <- NULL
    pespacers$revtranscript <- NULL
    pespacers$primer        <- NULL
    pespacers$pename <- pespacers$crisprname
    pespacers$type <- 'pespacer'

    nickingzone <- invertStrand(down_flank(pespacers, -3+40-5, -3+90+17))
    nickingspacers <- find_spacers(nickingzone, bsgenome, complement = FALSE)
    nickingspacers %<>% add_offtargets(bsgenome)
    nickingspacers %<>% filter_offtargets(groupby = 'pename')
    nickingspacers$type <- 'nickingspacer'
    
    # Plot
    if (plot){
        plotgr <- c(pespacers, nickingspacers)
        plotgr$type %<>% factor(c('pespacer', 'nickingspacer'))
        plot_intervals(plotgr, linetype_var = 'type', xref = 'pename', y = 'pename')
    }
    
    # Return
    nickingspacers
}
