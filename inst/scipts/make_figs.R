#=============
# TFBS
#=============

# Decide which to keep
require(multicrispr)
require(magrittr)
bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
bedfile  <- system.file('extdata/SRF.bed', package='multicrispr')
targets  <- bed_to_granges(bedfile, genome='mm10', plot = FALSE)
png('graphs/srf_karyogram.png')
plot_karyogram(targets, title = NULL)
dev.off()

extended <- extend(targets)
spacers  <- extended %>% 
            find_spacers(bsgenome) %>% 
            add_specificity(extended, bsgenome) %>% 
            add_efficiency(bsgenome, method = 'Doench2016', condaenv='azienv')
spacers %>% subset(seqnames == 'chr1') %>% 
            subset(specific == TRUE)   %>% 
            gr2dt() %>%
            extract(,.SD[.N>2] ,by = 'targetname') %>% 
            dt2gr(seqinfo(spacers))

targets %<>% subset(targetname %in% c('T0050', 'T0151'))
extended%<>% subset(targetname %in% c('T0050', 'T0151'))
spacers %<>% subset(targetname %in% c('T0050', 'T0151'))

# Create plots: original
p <- plot_intervals(targets) + 
    ggplot2::guides(color = FALSE, linetype = FALSE)
ggplot2::ggsave('graphs/srf.png', p, width = 2.5, height = 2)

targets$set <- 'target'
extended$set <- 'extended'
plotgr <- c(targets, extended)
plotgr$set %<>% factor(c('target', 'extended'))
p <- plot_intervals(plotgr, linetype_var = 'set') + 
    ggplot2::guides(color = FALSE, linetype = FALSE)
ggplot2::ggsave('graphs/srf_extended.png', p, width = 2.5, height = 2)


p <- plot_intervals(spacers) + 
    ggplot2::guides(color = FALSE)
ggplot2::ggsave('graphs/srf_spacers.png', p, width = 2.5, height = 2)

p <- plot_intervals(spacers, alpha_var = 'specific') + 
    ggplot2::scale_alpha_manual(values = c(`TRUE`=1, `FALSE`=0.20)) + 
    ggplot2::guides(color = FALSE, alpha = FALSE)
ggplot2::ggsave('graphs/srf_specific.png', p, width = 2.5, height = 2)


quantiles <- round(quantile(spacers$Doench2016, c(.33, .66, 1)), 2)
spacers$efficiency <- spacers$Doench2016 %>% cut(c(0, quantiles), labels = as.character(quantiles))
p <- plot_intervals(spacers, size_var = 'efficiency', alpha_var = 'specific') + 
    ggplot2::scale_size_manual(values = c(0.2, 1.5, 3) %>% set_names(quantiles)) + 
    ggplot2::scale_alpha_manual(values = c(`TRUE`=1, `FALSE`=0.25)) + 
    ggplot2::guides(color = FALSE, alpha = FALSE, size = FALSE)
ggplot2::ggsave('graphs/srf_efficient.png', p, width = 2.5, height = 2)

# PE
#=====
bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
gr <- char_to_granges(c(PRNP = 'chr20:4699600:+',             # snp
                        HBB  = 'chr11:5227002:-',              # snp
                        HEXA = 'chr15:72346580-72346583:-',    # del
                        CFTR = 'chr7:117559593-117559595:+'),  # ins
                       bsgenome)
p <- plot_intervals(subset(gr, targetname=='HBB'), facet_var = c('seqnames', 'targetname')) + 
    ggplot2::guides(color = FALSE)
ggplot2::ggsave('graphs/hbb.png', p, width=2.5, height=2)

extended <- extend_for_pe(gr, nrt = 48)
gr$set <- 'target'
grf <- grr <- gr
names(grf) %<>% paste0('_f')
names(grr) %<>% paste0('_r')
extended$set <- 'extension'
plotgr <- c(grf, grr, extended)
plotgr$set %<>% factor(c('target', 'extension'))

p <- plot_intervals(subset(plotgr, targetname=='HBB'), 
                    linetype_var = 'set', 
                    facet_var = c('seqnames', 'targetname')) + 
     ggplot2::guides(color = FALSE, linetype = FALSE)
ggplot2::ggsave('graphs/hbb_extended.png', p, width=2.5, height=2)


spacers <- gr %>% find_pe_spacers(bsgenome, nrt=48)

spacers %<>% add_genome_counts() %>%
            add_efficiency(bsgenome, 'Doench2016', condaenv = 'azienv')
spacers$specific <- spacers$G0==1
quantiles <- round(quantile(spacers$Doench2016, c(0.33, 0.66, 1)), 2)
spacers$efficiency <- cut(spacers$Doench2016, c(0, quantiles), labels = quantiles) %>% 
                        as.character()

spacers %<>% subset(targetname == 'HBB')
p <- plot_intervals(spacers, facet_var = c('seqnames', 'targetname')) + 
     ggplot2::guides(color = FALSE, size = FALSE)
ggplot2::ggsave('graphs/hbb_spacers.png', p, width=2.5, height=2)


p <- plot_intervals(spacers, facet_var = c('seqnames', 'targetname'), 
                    alpha_var = 'specific') + 
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.25)) + 
    ggplot2::guides(color = FALSE, size = FALSE, alpha = FALSE)
ggplot2::ggsave('graphs/hbb_specific.png', p, width=2.5, height=2)

p <- plot_intervals(spacers, facet_var = c('seqnames', 'targetname'), 
               alpha_var = 'specific', size_var  = 'efficiency') + 
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.25)) + 
    ggplot2::scale_size_manual(values = c(0.2, 1.5, 3) %>% set_names(quantiles)) + 
    ggplot2::guides(color = FALSE, size = FALSE, alpha = FALSE)
ggplot2::ggsave('graphs/hbb_efficient.png', p, width=2.5, height=2)


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
spacers <- find_pe_spacers(gr, bsgenome, nrt = 26)

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








