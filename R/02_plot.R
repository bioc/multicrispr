
#' Plot GRanges(List)
#' 
#' @details Visualize a GRanges(List) as karyogram (plot_karyogram), 
#' genome tracks (plot_tracks) or intervals (plot_intervals)
#' 
#' @param grangeslist \code{\link[GenomicRanges]{GRanges-class}} or 
#'                    \code{\link[GenomicRanges]{GRangesList-class}}
#' @param title plot title
#' @return list (plot_karyogram) or ggplot (plot_intervals)
#' @examples 
#' # Plot GRanges
#'   bedfile <-  system.file('extdata/SRF.bed',  package = 'multicrispr')
#'   granges   <- read_bed(bedfile, 'mm10' , plot = FALSE)
#'   plot_karyogram(granges)
#'   plot_intervals(granges)
#'   
#' # Plot GRangesList
#'   flanks  <- left_flank(granges)
#'   grangeslist <- GenomicRanges::GRangesList(sites = granges, flanks = flanks)
#'   plot_karyogram(grangeslist)
#'   plot_intervals(grangeslist)
#' @export
plot_karyogram <- function(grangeslist, title = unique(genome(grangeslist))){
    
    # Assert
    if (is(grangeslist, 'GRanges'))  grangeslist <- GRangesList(grangeslist)
    assert_is_identical_to_true(is(grangeslist, 'GRangesList'))
    
    # Extract
    genomeranges <- as(seqinfo(grangeslist)[seqlevelsInUse(grangeslist)], 
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


#' @rdname plot_karyogram
plot_tracks <- function(grangeslist){
    
    group <- . <-  NULL
    
    if (is(grangeslist, 'GRangesList')) granges <- unlist(grangeslist)
    genome  <- unique(genome(seqinfo(granges))); 
    assert_is_a_string(genome)
    chrom   <- unique(as.character(seqnames(granges)))[1]
    assert_is_a_string(chrom)
    
    # Find continuum groups
    granges$group <- findOverlaps(  granges, maxgap = 1, ignore.strand = TRUE, 
                                    select = 'first')
    
    # Plot
    coretracks <- list( ideogram = Gviz::IdeogramTrack(
                                        chromosome = chrom, 
                                        genome     = genome), 
                        genomeaxis = Gviz::GenomeAxisTrack())
    selectedgr <- subset(granges, group==1) %>% split(names(.))
    annottracks <- mapply(Gviz::AnnotationTrack, 
                          selectedgr, name = names(granges))
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
plot_intervals <- function(grangeslist){
    
    # Comply - Assert - Process
    contig <- group <- .N <- .SD <- tmp <- xstart <- xend <- y <- NULL
    assert_is_any_of(grangeslist, c('GRanges', 'GRangesList'))
    granges <-  if (is(grangeslist, 'GRangesList')){  unlist(grangeslist)
                } else {                              grangeslist    }

    # Find adjacent ranges    
    granges$contig <- findOverlaps(
        granges, maxgap = 1, select = 'first', ignore.strand = TRUE)
    granges %<>% extract(order(granges$contig))
    
    # Prepare plotdt
    plotdt <- data.table::as.data.table(granges)
    plotdt %>% extract( , 
            group := if (is.null(names(granges))) 'ranges' else names(granges))
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
        ggplot2::theme_bw()

    # Print and return
    print(p)
    p
    
}

