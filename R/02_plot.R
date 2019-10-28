
#' Karyo/Interval Plot GRanges(List)
#' @param grangeslist \code{\link[GenomicRanges]{GRanges-class}} or 
#'                    \code{\link[GenomicRanges]{GRangesList-class}}
#' @param title plot title
#' @return list (plot_karyogram) or ggplot (plot_intervals)
#' @seealso  \code{\link[karyoploteR]{plotKaryotype}}, around which this 
#'           function wraps
#' @examples 
#' # Plot GRanges
#'   bedfile <-  system.file('extdata/SRF.bed',  package = 'multicrispr')
#'   bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'   gr <- bed_to_granges(bedfile, bsgenome, plot = FALSE)
#'   plot_karyogram(gr)
#'   plot_intervals(gr)
#'   
#' # Plot GRangesList
#'   flanks  <- left_flank(gr)
#'   grangeslist <- GenomicRanges::GRangesList(sites = gr, flanks = flanks)
#'   plot_karyogram(grangeslist)
#'   plot_intervals(grangeslist)
#' @export
plot_karyogram <- function(
    grangeslist, 
    title = unique(GenomeInfoDb::genome(grangeslist))
){
    
    # Assert
    . <- NULL
    if (is(grangeslist, 'GRanges'))  grangeslist <- GRangesList(grangeslist)
    assertive.types::assert_is_all_of(grangeslist, 'GRangesList')
    
    # Extract relevant chromosomes and order them
    chroms <- union(GenomeInfoDb::seqlevelsInUse(grangeslist), 
                    canonicalseqlevels(grangeslist))
    stri_extract <- function(stri, pattern){
        stri %>% extract(stringi::stri_detect_regex(., pattern)) 
    }
    chrs1  <- chroms %>% stri_extract('^(chr)?[0-9]$')    %>% sort()
    chrs2  <- chroms %>% stri_extract('^(chr)?[0-9]{2}$') %>% sort()
    chrsXY <- chroms %>% stri_extract('^(chr)?[XY]$')     %>% sort()
    chrsM  <- chroms %>% stri_extract('^(chr)?MT?$')
    orderedchrs <-  c(chrs1, chrs2, chrsXY, chrsM) %>% 
                    c(sort(setdiff(chroms, .)))
    genomeranges <- as(GenomeInfoDb::seqinfo(grangeslist)[orderedchrs],
                        "GRanges")

    # Color
    n <- length(grangeslist)
    if (n>0){
        hues <- seq(15, 375, length = n + 1)
        colors  <-  grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)]
    }
    
    # Plot
    kp <- karyoploteR::plotKaryotype(genomeranges, main = title)
    for (i in seq_len(n)){
        karyoploteR::kpPlotRegions(kp, grangeslist[[i]], col = colors[i])
    }
    
    # Add legend
    if (assertive.properties::has_names(grangeslist)){
        graphics::legend('right', fill = colors, legend = names(grangeslist))
    }

}


plot_tracks <- function(grangeslist){
    
    group <- . <-  NULL
    
    if (is(grangeslist, 'GRangesList')) gr <- unlist(grangeslist)
    genome  <- unique(genome(seqinfo(gr)))
    assertive.types::assert_is_a_string(genome)
    chrom   <- unique(as.character(seqnames(gr)))[1]
    assertive.types::assert_is_a_string(chrom)
    
    # Find continuum groups
    gr$group <- findOverlaps(  gr, maxgap = 1, ignore.strand = TRUE, 
                                    select = 'first')
    
    # Plot
    coretracks <- list( ideogram = Gviz::IdeogramTrack(
                                        chromosome = chrom, 
                                        genome     = genome), 
                        genomeaxis = Gviz::GenomeAxisTrack())
    selectedgr   <- subset(gr, group==1) %>% split(names(.))
    annottracks  <- mapply( Gviz::AnnotationTrack, 
                            selectedgr, name = names(gr))
    Gviz::plotTracks(c(coretracks, annottracks), 
                    background.title = 'gray40', 
                    add = TRUE)

}


to_megabase <- function(y){
    z <- vector('character', length(y))
    
    i <- y>1e6
    z[i] <- paste0(round(y[i]*1e-6), ' M')
    
    i <- y>1e3 & y<=1e6
    z[i] <- paste0(round(y[i]*1e-3), ' K')
    
    i <- y<=1e3
    z[i] <- paste0(round(y[i]), 'b')
    z %>% set_names(names(y))
}


#' @rdname plot_karyogram
#' @export
plot_intervals <- function(grangeslist, title = NULL){
    
    # Comply - Assert - Process
    contig <- group <- .N <- .SD <- tmp <- xstart <- xend <- y <- NULL
    assertive.types::assert_is_any_of(grangeslist, c('GRanges', 'GRangesList'))
    gr <-  if (is(grangeslist, 'GRangesList')){  unlist(grangeslist)
                } else {                              grangeslist    }

    # Find adjacent ranges    
    gr$contig <- findOverlaps(
        gr, maxgap = 1, select = 'first', ignore.strand = TRUE)
    gr %<>% extract(order(gr$contig))
    
    # Prepare plotdt
    plotdt <- data.table::as.data.table(gr)
    plotdt %>% extract( , 
            group := if (is.null(names(gr))) 'ranges' else names(gr))
    plotdt <- plotdt[ , .SD[contig %in% c(min(contig), max(contig)) ], 
                        by = c('seqnames')]
    plotdt %<>% extract(order(seqnames, start))
    plotdt %>%  extract(, y      := min(start), by = 'contig')
    plotdt %>%  extract(, y      := factor(round(y*1e-6)))
    plotdt %>%  extract(, xstart := start-min(start), by = 'contig')
    plotdt %>%  extract(, xend   := xstart + width)
    plotdt %>%  extract(strand=='-', tmp    := xend)
    plotdt %>%  extract(strand=='-', xend   := xstart)
    plotdt %>%  extract(strand=='-', xstart := tmp)
    plotdt %>%  extract(, tmp := NULL)
    
    # Plot
    p <-ggplot2::ggplot(plotdt) + 
        ggplot2::facet_wrap(~ seqnames, scales = 'free') + 
        ggplot2::geom_segment(
            ggplot2::aes(x = xstart, xend = xend, y = y, yend = y, color=group),
            arrow = grid::arrow(length = ggplot2::unit(0.1, "inches"))) + 
        ggplot2::theme_bw() + 
        ggplot2::xlab('Bases') + 
        ggplot2::ylab('Megabases') + 
        ggplot2::ggtitle(title)

    # Print and return
    print(p)
    p
    
}

