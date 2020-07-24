require(magrittr)
require(multicrispr)

#=============
# TFBS
#=============

# Decide which to keep
reticulate::use_condaenv('azienv')
bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
bedfile  <- system.file('extdata/SRF.bed', package='multicrispr')
targets  <- multicrispr::bed_to_granges(bedfile, genome='mm10', plot = FALSE)
#png('graphs/srf_karyogram.png')
#multicrispr::plot_karyogram(targets, title = NULL)
#dev.off()

extended <- extend(targets, -22, +22)

# targets['T0151'] %>% plot_intervals()
# targets['T0151'] %>% extend()     %>% plot_intervals()
# targets['T0151'] %>% up_flank()   %>% plot_intervals()
# targets['T0151'] %>% down_flank() %>% plot_intervals()

spacers  <- extended %>% find_spacers(bsgenome)
spacers %<>% add_offtargets(bsgenome, extended)
spacers %<>% add_efficiency(bsgenome, method = 'Doench2016')

spacers %>% subset(seqnames == 'chr1') %>% 
            subset(specific == TRUE)   %>% 
            gr2dt() %>%
            extract(,.SD[.N>2] ,by = 'targetname') %>% 
            dt2gr(seqinfo(spacers))

spacers2 <- spacers
spacers2 %<>% subset(specific==TRUE)
spacers2$score <- spacers2$name <- spacers2$crisprcontext <- spacers2$specific <- NULL
spacers2

targets %<>% subset(targetname %in% c('T0050')) # T0151
extended%<>% subset(targetname %in% c('T0050'))
spacers %<>% subset(targetname %in% c('T0050'))
spacers %>% extract(c(length(.), 1))

# Create plots: original
p <- plot_intervals(targets, color_var = 'targetname') +
    ggplot2::scale_color_manual(values = c(T0050 = "#00BFC4")) +
    ggplot2::guides(color = FALSE, linetype = FALSE)
ggplot2::ggsave('graphs/srf01.pdf',   p, width =2.2, height = 1.5, device = grDevices::cairo_pdf)
#ggplot2::ggsave('graphs/srf01.png', p, width = 2.2, height = 1.5)


p <- plot_intervals(extended) + 
    ggplot2::scale_color_manual(values = c(T0050 = "#00BFC4")) +
    ggplot2::guides(color = FALSE, linetype = FALSE)
ggplot2::ggsave('graphs/srf02_extended.pdf', p, width =2.2, height = 1.5, device = grDevices::cairo_pdf)
#ggplot2::ggsave('graphs/srf02_extended.png', p, width = 2.2, height = 1.5)

p <- plot_intervals(spacers) + 
    ggplot2::scale_color_manual(values = c(T0050 = "#00BFC4")) +
    ggplot2::guides(color = FALSE)
ggplot2::ggsave('graphs/srf03_spacers.pdf', p, width =2.2, height = 1.5, device = grDevices::cairo_pdf)
#ggplot2::ggsave('graphs/srf03_spacers.png', p, width = 2.2, height = 1.5)

p <- plot_intervals(spacers, alpha_var = 'specific') + 
    ggplot2::scale_color_manual(values = c(T0050 = "#00BFC4")) +
    ggplot2::scale_alpha_manual(values = c(`TRUE`=1, `FALSE`=0.20)) + 
    ggplot2::guides(color = FALSE, alpha = FALSE)
ggplot2::ggsave('graphs/srf04_specific.pdf', p, width = 2.2, height = 1.5, device = grDevices::cairo_pdf)
#ggplot2::ggsave('graphs/srf04_specific.png', p, width = 2.2, height = 1.5)

quantiles <- round(quantile(spacers2$Doench2016, c(.33, .66, 1)), 2)
spacers$efficiency <- spacers$Doench2016 %>% cut(c(0, quantiles), labels = as.character(quantiles))
p <- plot_intervals(spacers, size_var = 'efficiency', alpha_var = 'specific') + 
    ggplot2::scale_color_manual(values = c(T0050 = "#00BFC4")) +
    ggplot2::scale_size_manual(values = c(0.2, 1.5, 3) %>% set_names(quantiles)) + 
    ggplot2::scale_alpha_manual(values = c(`TRUE`=1, `FALSE`=0.25)) + 
    ggplot2::guides(color = FALSE, alpha = FALSE, size = FALSE)
ggplot2::ggsave('graphs/srf05_efficient.pdf', p, width = 2.2, height = 1.5, device = grDevices::cairo_pdf)
#ggplot2::ggsave('graphs/srf05_efficient.png', p, width = 2.2, height = 1.5)

# PE
#=====
require(magrittr)
require(multicrispr)
require(ggplot2)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
gr <- char_to_granges(c(PRNP='chr20:4699600:+'), bsgenome)
find_primespacers(gr, bsgenome, edits = 'T')
plot_intervals(gr)
blanken <- function(
    p, 
    axis.text.x        = element_blank(),
    axis.text.y        = element_blank(),
    axis.ticks.x       = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.x = element_line(size=.5),
    panel.grid.minor.x = element_line(size=.5)
    
){
    p + 
    theme(  strip.background   = element_blank(), 
            strip.text.x       = element_blank(), 
            axis.text.x        = axis.text.x,
            axis.text.y        = axis.text.y,
            axis.ticks.x       = axis.ticks.x,
            axis.ticks.y       = axis.ticks.y,
            panel.grid.major.y = element_blank(), 
            panel.grid.minor.y = element_blank(), 
            panel.grid.major.x = panel.grid.major.x, 
            panel.grid.minor.x = panel.grid.minor.x, 
            plot.background    = element_rect(fill = 'transparent', colour=NA)) + 
    guides(color = FALSE, size = FALSE, linetype = FALSE, alpha = FALSE)
}

extended <- extend_for_pe(gr)
plot_intervals(extended) %>% blanken()
ggsave('../graphs/prnp02_extended.pdf', width=1.3, height=0.6, device = grDevices::cairo_pdf, bg = 'transparent')

spacers <- find_primespacers(gr, bsgenome, ontargets = 'Doench2016')
(plot_intervals(spacers, alpha_var = 'type', size_var = NULL) + 
scale_alpha_manual(values = c(`spacer` = 1, `3 extension` = 1, `nicking spacer` = 0))) %>% blanken()
ggplot2::ggsave('../graphs/prnp03_primespacers.pdf', width=1.3, height=0.9, device = grDevices::cairo_pdf, bg = 'transparent')

plot_intervals(spacers, alpha_var = NULL, size_var = NULL) %>% blanken(axis.text.y = element_text())
ggplot2::ggsave('../graphs/prnp04_nickspacers.pdf',  width=1.85, height=0.9, device = grDevices::cairo_pdf, bg = 'transparent')

plot_intervals(spacers) %>% blanken(axis.text.x = element_text())
ggplot2::ggsave('../graphs/prnp05_onofftargets.pdf', width=1.3, height=1.0, device = grDevices::cairo_pdf, bg = 'transparent')



spacers %<>% add_ontargets(bsgenome, 'Doench2016')
plot_intervals(spacers) %>% blanken()

quantiles <- round(quantile(spacers$Doench2016, c(0.33, 0.66, 1)), 2)
spacers$efficiency <- cut(spacers$Doench2016, c(0, quantiles), labels = quantiles) %>% 
                        as.character()
p <- plot_intervals(spacers, facet_var = 'seqnames') + 
    ggplot2::guides(color = FALSE, size = FALSE, alpha = FALSE)
ggplot2::ggsave('graphs/hbb05_efficient.pdf', p, width=2.2, height=1.8, device = grDevices::cairo_pdf)


gr %<>% extract('HBB')
plot_intervals(gr, facet_var = c('targetname', 'seqnames'), color_var = NULL) + 
    #ggplot2::guides(color=FALSE) + 
    ggplot2::xlab(NULL) + ggplot2::ylab(NULL) + 
    ggplot2::scale_x_continuous(name = NULL, breaks = c(0, 1)) + 
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank())


# Extend
find_spacers(extend_for_pe(gr, bsgenome), bsgenome, complement = FALSE)

grplot <- c(gr, grext)
grplot$set <- c('target', '+', '-')
plot_intervals(grplot, facet_var = c('targetname', 'seqnames'), yby = 'set', color_var = NULL, linetype = 'set') + 
    ggplot2::xlab(NULL) + ggplot2::ylab(NULL)  + 
    #ggplot2::guides(color=FALSE) + 
    ggplot2::scale_linetype_manual(values = c(target = 'solid', `+` = 'longdash', `-` = 'longdash')) + 
    #ggplot2::scale_color_manual(values = c(target = '#F8766D', `+` = 'black' , `-` = 'black')) + 
    ggplot2::guides(linetype=FALSE) + 
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                   panel.grid.minor = ggplot2::element_blank())
grDevices::hcl(h = seq(15, 375, length = n + 1), l = 65, c = 100)[1:3]

    #ggplot2::theme( legend.position   = c(0.25,0.73), 
    #                legend.title      = ggplot2::element_blank(), 
    #                legend.key        = ggplot2::element_blank(), 
    #                legend.background = ggplot2::element_blank())

# Find pe sites
spacers <- find_primespacers(gr, bsgenome, nrt = 26)

# Filter for specificity
genomedir <- '~/.multicrispr/bowtie/genome/BSgenome.Hsapiens.UCSC.hg38'
outdir <- '~/.multicrispr/bowtie'
filter_prime_specific(spacers, genomedir, outdir)



# quantiles <- round(quantile(spacers$Doench2016, c(.33, .66, 1)), 2)
# spacers$efficiency <- spacers$Doench2016 %>% cut(c(0, quantiles), labels = as.character(quantiles))
# p <- plot_intervals(spacers, size_var = 'efficiency') + 
#     ggplot2::scale_size_manual(values = c(0.2, 1.5, 4) %>% set_names(quantiles)) + 
#     ggplot2::guides(color = FALSE, alpha = FALSE, size = FALSE) 
# ggplot2::ggsave('graphs/srf_efficient.png', p, width = 3, height = 2.5)

p <- plot_intervals(spacers, size_var = 'efficiency', alpha_var = 'specific') + 
    ggplot2::scale_size_manual(values = c(0.2, 1.5, 3) %>% set_names(quantiles)) + 
    ggplot2::scale_alpha_manual(values = c(`TRUE`=1, `FALSE`=0.25)) + 
    ggplot2::guides(color = FALSE, alpha = FALSE, size = FALSE)
ggplot2::ggsave('graphs/srf_efficient.png', p, width = 2.5, height = 2)


# Process
png('graphs/srf_extended.png')
targets %<>% extend(plot = TRUE)
dev.off()


targets %>% find_spacers(bsgenome) %>% filter_target_specific(targets, bsgenome) %>% filter_efficient()

targets %<>% subset(targetname %in% c('T0018', 'T1042'))
plot_intervals(targets, yby = 'targetname')

# Extend
grext <- extend(gr, plot = TRUE)
gr$set <- 'target'
grext$set <- 'extension'
grplot <- c(gr, grext)
grplot$set %<>% factor(c('target', 'extension'))
plot_intervals(grplot, color_var = 'targetname', linetype = 'set')

# Find spacers
spacers <- find_spacers(grext, bsgenome)
spacers %>% subset(targetname %in% c('T0018', 'T0151')) %>% plot_intervals(yby = 'crisprname', color_var = 'targetname')

# Filter for specificity
tdir <- index_targets(grext, bsgenome)
gdir <- '~/.multicrispr/bowtie/genome/BSgenome.Mmusculus.UCSC.mm10'
outdir <- '~/.multicrispr/bowtie'
specificspacers <- filter_target_specific(spacers, tdir, gdir, outdir)
specificspacers %>% subset(targetname %in% c('T0018', 'T0151')) %>% plot_intervals(yby = 'crisprname', color_var = 'targetname')
    
# Filter for 
specificspacers %<>% score_spacers(bsgenome, 'Doench2016', condaenv = 'azimuthenv')
specificspacers %>% subset(targetname %in% c('T0018', 'T0151') & Doench2016 > 0.4)

    plot_intervals(specificspacers, yby = 'crisprname')





targets  <- extend(bed_to_granges(bedfile, genome='mm10', plot = FALSE))
bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
spacers <- find_spacers(targets, bsgenome)
indexedgenome <- '~/.multicrispr/bowtie/BSgenome.Mmusculus.UCSC.mm10'
specific <- filter_target_specific(spacers, targets, bsgenome, indexedgenome)
specific %<>% score_spacers(bsgenome, 'Doench2016', condaenv = 'azimuthenv')
efficient <- specific %<>% subset(Doench2016>0.4)


targets$set   <- 'A'
spacers$set   <- 'B'
specific$set  <- 'C'
efficient$set <- 'D'

plotgr <- c(targets, spacers, specific, efficient)

    plotgr %<>% sort(ignore.strand = TRUE)
    #plotgr %<>% subset(seqnames == 'chr1')
    plotdt <- data.table::as.data.table(plotgr)
    plotdt[ , bin := cut(start, 10, labels = FALSE) ]
    plotdt[ , y := scales::number(min(start), big.mark = ','), by = 'bin']
    plotdt$y %<>% factor(unique(.))
    
    plotdt[strand=='+', x    := start]
    plotdt[strand=='+', xend := end]
    plotdt[strand=='-', x    := end]
    plotdt[strand=='-', xend := start]
    
    # plotdt[strand=='+', x    := start - min(start), by = 'y' ]
    # plotdt[strand=='+', xend := end   - min(start), by = 'y' ]
    # plotdt[strand=='-', x    := end   - min(start), by = 'y' ]
    # plotdt[strand=='-', xend := start - min(start), by = 'y' ]
    
    require(ggplot2)
    ggplot(plotdt, aes_string(x = 'x', xend = 'xend', y = 'seqnames', yend = 'seqnames', color = 'set')) + 
    geom_segment(arrow = arrow(length = unit(0.1, "inches"))) + 
    scale_x_continuous(labels = scales::comma) + 
    xlab(NULL) + ylab(NULL) + theme_bw() + 
    scale_color_manual(values=c(A='lavenderblush', B='gray90', C='gray50', D='green4'))
    #guides(color = FALSE)








