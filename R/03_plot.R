
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
#'     flanks  <- up_flank(gr, stranded=FALSE)
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


strsplitextract <- function(x, split, fixed = FALSE, i){
    vapply(strsplit(x, split = split, fixed=fixed), extract, character(1), i)
}



plot_intervals_engine <- function(
    gr, xref, y, nperchrom, nchrom , color_var, facet_var, linetype_var, 
    size_var, alpha_var, title, scales
){
# Comply
    edge <- targetname <- NULL
# Assert, Import, Comply
    assert_is_all_of(gr, 'GRanges')
    if (assertive::is_empty(gr)) return(invisible(NULL))
    if (!is.null(color_var)) assert_is_a_string(color_var)
    assert_is_subset(color_var, names(as.data.table(gr)))
    contig <- .N <- .SD <- seqnames <- start <- NULL
    strand <- tmp <- width <- xstart <- xend <- . <- NULL
# Prepare plotdt
    plotdt <- prepare_plot_intervals(
                gr, xref, y, nperchrom, nchrom, alpha_var, size_var)
# Core Ranges
    p <-ggplot( plotdt, 
                aes_string(x = 'xstart', xend = 'xend', y = 'y', yend = 'y', 
                            color = color_var, linetype = linetype_var, 
                            size = size_var, alpha = alpha_var)) + 
        facet_wrap(facet_var, scales = scales) + 
        geom_segment(arrow = arrow(length = unit(0.1, "inches")))
    if (!is.null(alpha_var)) p <- p + 
        scale_alpha_manual(values = c(`0` = 1, `1+` = 0.3))
    if (!is.null(size_var))  p <- p +
        scale_size_manual(values = c(`0+` = 0.1, `0.3+` = 1, `0.5+` = 2))
# Targets
    if (all(c('targetstart', 'targetend') %in% names(mcols(gr)))){
        p <-p + geom_point(aes_string(
                    x = 'xtargetstart', y = 'y'), shape = '|', size = 4) +
                geom_point(aes_string(
                    x = 'xtargetend',   y = 'y'), shape = '|', size = 4)}
# Extensions
    p <- p + theme_bw()  +  xlab(NULL)  +  ylab(NULL)  +  ggtitle(title)
# Return
    p # print(p)
}


#' Interval plot GRanges
#' 
#' @param gr          \code{\link[GenomicRanges]{GRanges-class}}
#' @param xref        gr var used for scaling x axis
#' @param y           'names' (default) or name of gr variable
#' @param nperchrom    number (default 1): n head (and n tail) targets 
#'                     shown per chromosome
#' @param nchrom       number (default 6) of chromosomes shown
#' @param color_var   'seqnames' (default) or other gr variable
#' @param linetype_var NULL (default) or gr variable mapped to linetype
#' @param size_var     NULL (default) or gr variable mapped to size
#' @param facet_var    NULL(default)  or gr variable mapped to facet
#' @param alpha_var    NULL or gr variable mapped to alpha
#' @param title        NULL or string: plot title
#' @param scales       'free', 'fixed', etc
#' @return ggplot object
#' @seealso  \code{\link{plot_karyogram}}
#' @examples 
#' # SRF sites
#'     require(magrittr)
#'     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
#'     bedfile <-  system.file('extdata/SRF.bed',  package = 'multicrispr')
#'     targets   <- bed_to_granges(bedfile, 'mm10', plot = FALSE)
#'     plot_intervals(targets)
#'     
#' # PE targets
#'     bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
#'     gr <- char_to_granges(c(PRNP = 'chr20:4699600:+',
#'                             HBB  = 'chr11:5227002:-',
#'                             HEXA = 'chr15:72346580-72346583:-',
#'                             CFTR = 'chr7:117559593-117559595:+'), 
#'                           bsgenome)
#'     spacers <- find_pe_spacers(gr, bsgenome, plot = FALSE)
#'     plot_intervals(gr)
#'     plot_intervals(extend_for_pe(gr))
#'     plot_intervals(spacers)
#'     
#' # Empty gr
#'     plot_intervals(GenomicRanges::GRanges())
#' @export
plot_intervals <- function(
    gr, xref = 'targetname', y = default_y(gr), nperchrom = 2, nchrom = 4, 
    color_var = 'targetname', facet_var = 'seqnames', 
    linetype_var = default_linetype(gr), size_var = default_size_var(gr), 
    alpha_var = default_alpha_var(gr), title = NULL, scales= 'free'
){
    gr$type <- 'spacer'
    gr$type %<>% factor(c("spacer", "3' extension", "nicking spacer"))
    seqinfo1 <- seqinfo(gr)
    plotgr <- gr

    if ('crisprextension' %in% names(mcols(gr))){
        extgr <- GRanges(gr$crisprextrange, seqinfo=seqinfo1)
        mcols(extgr) <- mcols(gr)
        extgr$type <- "3' extension"
        plotgr %<>% c(extgr)
    }
    
    if ('nickrange' %in% names(mcols(gr))){
        nickdt <- gr2dt(gr)
        nickdt %<>% extract(complete.cases(nickrange))
        nickdt %<>%  tidyr::separate_rows(tidyselect::starts_with('nick'), sep = ';')
        nickgr <- GRanges(nickdt$nickrange, seqinfo = seqinfo1)
        nickgr$off <- as.numeric(nickdt$nickoff)
        mcols(nickgr) <- mcols(dt2gr(nickdt, seqinfo = seqinfo1))
        nickgr$type <- 'nicking spacer'
        plotgr %<>% c(nickgr)
    }
    
    plot_intervals_engine(  plotgr, 
                            xref         = xref, 
                            y            = y, 
                            nperchrom    = nperchrom,
                            nchrom       = nchrom,
                            color_var    = color_var, 
                            facet_var    = facet_var,
                            linetype     = linetype_var, 
                            size_var     = size_var,
                            alpha_var    = alpha_var,
                            title        = title, 
                            scales       = scales)
    
}

default_linetype <- function(gr){
    if (any(c('crisprextension', 'nickrange') %in% names(mcols(gr)))){
        'type'
    } else {
        NULL
    }
}

default_y <- function(gr){
    if ('crisprname' %in% names(mcols(gr)))  'crisprname' else 'names'
    
}

default_alpha_var <- function(gr){
    if ('off' %in% names(mcols(gr))) 'off' else NULL
}

default_size_var <- function(gr){
    if ('Doench2016' %in% names(mcols(gr))) 'Doench2016' else NULL
}


head_tail <- function(x, n){
    idx <- x %in% unique( c(head(x, ceiling(n/2)), tail(x, floor(n/2))))
    x[idx]
}

prepare_plot_intervals <- function(
    gr, xref, y, nperchrom, nchrom, alpha_var, size_var
){
# Comply
    edge <- targetname <- xstart <- xend <- width <- NULL
    targetstart <- targetend <- xtargetstart <- xtargetend <- NULL
    extstart <- crisprprimer <- crisprtranscript <- crisprextension <- tmp <- NULL
# Prepare data.table. Select chromosomes/targets to plot.
    plotdt <- data.table::as.data.table(gr) %>% cbind(names = names(gr))
    plotdt %<>% extract(order(seqnames, start))
    plotdt$seqnames %<>% droplevels()
    headtailchroms <- head_tail(levels(plotdt$seqnames), nchrom)
    plotdt %<>% extract(headtailchroms, on = 'seqnames')
    plotdt$seqnames %<>% factor(headtailchroms)
    plotdt %<>% extract( # targets
        , .SD[targetname %in% head_tail(unique(targetname), nperchrom)],
        by = 'seqnames')
# Main ranges
    plotdt %>%  extract(, y      := min(start), by = y)
    plotdt %>%  extract(, y      := factor(format(y, big.mark = " ")))
    plotdt %>%  extract(, xstart := start-min(start), by = xref)
    plotdt %>%  extract(, xend   := xstart + width)
# Target marks
    if (all(c('targetstart', 'targetend') %in% names(mcols(gr)))){
        plotdt %>% extract(, xtargetstart := xstart + targetstart-start)
        plotdt %>% extract(, xtargetend   := xend   + targetend-end  )}
# Flip for arrow direction    
    plotdt[strand=='-', tmp    := xend]
    plotdt[strand=='-', xend   := xstart]
    plotdt[strand=='-', xstart := tmp]
    plotdt[           , tmp      := NULL]
# Alpha and Size
    if (!is.null(alpha_var)) plotdt[[alpha_var]] %<>% cut(
                                    c(-Inf, 0, Inf), c('0', '1+'))
    if (!is.null(size_var))  plotdt[[size_var]]  %<>% cut(
                                    c(-Inf, 0.3, 0.5, Inf), 
                                c('0+', '0.3+', '0.5+'))
# Return
    plotdt
}
