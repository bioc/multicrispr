# Load
require(magrittr)
require(multicrispr)
devtools::load_all()

# PE
#=====
bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
gr <- char_to_granges(c(PRNP = 'chr20:4699600:+',             # snp
                        HBB  = 'chr11:5227002:-',              # snp
                        HEXA = 'chr15:72346580-72346583:-',    # del
                        CFTR = 'chr7:117559593-117559595:+'),  # ins
                       bsgenome)
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


#=============
# TFBS
#=============

# Read
bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
targets  <- bed_to_granges(bedfile, genome='mm10', plot = FALSE)
png('graphs/srf.png')
plot_karyogram(targets, title = NULL)
dev.off()

# Process
targets %<>% extend()
bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
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








