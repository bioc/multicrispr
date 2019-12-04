
#' Karyo/Interval Plot GRanges(List)
#' @param grlist \code{\link[GenomicRanges]{GRanges-class}}
#' @param title plot title
#' @return list
#' @seealso  \code{\link{plot_intervals}}
#' @examples 
#' # Plot GRanges
#'     bedfile <-  system.file('extdata/SRF.bed',  package = 'multicrispr')
#'     gr <- bed_to_granges(bedfile, 'mm10', plot = FALSE)
#'     plot_karyogram(gr)
#'   
#' # Plot GRangesList
#'     flanks  <- left_flank(gr)
#'     grlist <- GenomicRanges::GRangesList(sites = gr, flanks = flanks)
#'     plot_karyogram(grlist)
#' @export
plot_karyogram <- function(
    grlist, 
    title = unique(genome(grlist))
){
    
    # Assert
    . <- NULL
    if (methods::is(grlist, 'GRanges')){
        grlist <- GenomicRanges::GRangesList(grlist)
    }
    assert_is_all_of(grlist, 'GRangesList')
    
    # Extract relevant chromosomes and order them
    chroms <- union(seqlevelsInUse(grlist), standardChromosomes(grlist))
    stri_extract <- function(stri, pattern){
        stri %>% extract(stri_detect_regex(., pattern)) 
    }
    chrs1  <- chroms %>% stri_extract('^(chr)?[0-9]$')    %>% sort()
    chrs2  <- chroms %>% stri_extract('^(chr)?[0-9]{2}$') %>% sort()
    chrsXY <- chroms %>% stri_extract('^(chr)?[XY]$')     %>% sort()
    chrsM  <- chroms %>% stri_extract('^(chr)?MT?$')
    orderedchrs <-  c(chrs1, chrs2, chrsXY, chrsM) %>% 
                    c(sort(setdiff(chroms, .)))
    genomeranges <- as(seqinfo(grlist)[orderedchrs], "GRanges")

    # Color
    n <- length(grlist)
    if (n>0){
        hues <- seq(15, 375, length = n + 1)
        colors  <-  grDevices::hcl(h = hues, l = 65, c = 100)[seq_len(n)]
    }
    
    # Plot
    kp <- karyoploteR::plotKaryotype(genomeranges, main = title)
    for (i in seq_len(n)){
        karyoploteR::kpPlotRegions(kp, grlist[[i]], col = colors[i])
    }
    
    # Add legend
    if (has_names(grlist)){
        graphics::legend('right', fill = colors, legend = names(grlist))
    }

}


plot_tracks <- function(grlist){
    
    group <- . <-  NULL
    
    if (methods::is(grlist, 'GRangesList')) gr <- unlist(grlist)
    genome  <- unique(genome(seqinfo(gr)))
    assert_is_a_string(genome)
    chrom   <- unique(as.character(seqnames(gr)))[1]
    assert_is_a_string(chrom)
    
    # Find continuum groups
    gr$group <- GenomicRanges::findOverlaps(
                    gr, maxgap = 1, ignore.strand = TRUE, select = 'first')
    
    # Plot
    coretracks <- list( ideogram = Gviz::IdeogramTrack(
                                        chromosome = chrom, 
                                        genome     = genome), 
                        genomeaxis = Gviz::GenomeAxisTrack())
    selectedgr   <- subset(gr, group==1) %>% split(names(.))
    annottracks  <- mapply( Gviz::AnnotationTrack, selectedgr, name = names(gr))
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


#' Interval plot GRanges
#' @param gr        \code{\link[GenomicRanges]{GRanges-class}}
#' @param color_var string: mcol mapped to color.
#' @param title     plot title
#' @return ggplot object
#' @seealso  \code{\link{plot_karyogram}}
#' @examples 
#' # SRF sites
#'     require(magrittr)
#'     bedfile <-  system.file('extdata/SRF.bed',  package = 'multicrispr')
#'     sites   <- bed_to_granges(bedfile, 'mm10', plot = FALSE)
#'     plot_intervals(sites)
#'     
#'     flanks  <- left_flank(sites)
#'     sites$color <- 'sites'
#'     flanks$color <- 'flanks'
#'     plot_intervals(c(sites, flanks))
#'     
#' # PE targets
#'     bs <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#'     sites <- GenomicRanges::GRanges(
#'                seqnames= c(PRNP = 'chr20:4699600',
#'                            HBB  = 'chr11:5227002',
#'                            HEXA = 'chr15:72346580-72346583',
#'                            CFTR = 'chr7:117559593-117559595'),
#'                strand   = c(PRNP = '+', HBB = '-', HEXA = '-', CFTR = '+'), 
#'                seqinfo  = BSgenome::seqinfo(bs))
#'     sites$color <- names(sites)
#'     plot_intervals(sites)
#' @export
plot_intervals <- function(
    gr, 
    color_var = if ('color' %in% names(mcols(gr))) 'color' else 'seqnames', 
    title = NULL
){
    
    # Assert, Import, Comply
    assert_is_all_of(gr, 'GRanges')
    assert_is_subset(color_var, names(as.data.table(gr)))
    contig <- group <- .N <- .SD <- seqnames <- start <- NULL
    strand <- tmp <- width <- xstart <- xend <- y <- NULL

    # Find adjacent ranges    
    gr$contig <- GenomicRanges::findOverlaps(
                    gr, maxgap = 1, select = 'first', ignore.strand = TRUE)
    gr %<>% extract(order(gr$contig))
    
    # Prepare plotdt
    plotdt <- as.data.table(gr)
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
    p <-ggplot( plotdt, 
                aes_string( x = 'xstart', xend = 'xend', y = 'y', yend = 'y', 
                            color = color_var)) + 
        facet_wrap(~ seqnames, scales = 'free') + 
        geom_segment(arrow = arrow(length = unit(0.1, "inches"))) + 
        theme_bw() + 
        xlab('Bases') + ylab('Megabases') + ggtitle(title)

    # Print and return
    print(p)
    p
    
}
