# Load functions
    require(reticulate)
    use_condaenv('azienv')
    import('azimuth')
    # require(AnnotationDbi)
    # require(BSgenome)
    # require(multicrispr)
    # require(magrittr)
    # require(biomaRt)
    # require(data.table)
    # require(stringi)

# Reconstruct Brunello spacers . Transform into exons
    brunello <- multicrispr:::reconstruct_brunello_spacers()
    brunello_spacers <- multicrispr:::validate_brunello_spacers(brunello)
    brunello_exons <- multicrispr:::find_enclosing_exons(brunello)
    brunello %<>% extract(brunello_exons$targetname)

# Validate
    bsgenome <- BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
    exonspacers <- find_spacers(brunello_exons, bsgenome, plot = FALSE) # 2 583 378
    exonspacers %<>% add_efficiency(bsgenome, 'Doench2016', chunksize=100000, plot = FALSE)
    bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    seqlevelsStyle(exonspacers) <- 'UCSC'
    exonspacers %<>% add_genome_counts(
                        bsgenome, outdir = '../multicrisprout/brunello', 
                        indexedgenomesdir = '../multicrisprout/indexedgenomes', 
                        mismatches = 1)
    seqlevelsStyle(exonspacers) <- 'NCBI'
    #saveRDS(exonspacers, '../multicrisprout/brunello/exonspacers.rds')
    #exonspacers <- readRDS('../multicrisprout/brunello/exonspacers.rds')
    
# Plot    
    # Venn
        brunellocoords   <- unname(as.character(granges(brunello)))
        exonspacercoords <- unname(as.character(granges(exonspacers)))
        x <- list(multicrispr = exonspacercoords, Brunello = brunellocoords)
        p_brunello_venn <- autonomics::plot_venn(x, title = 'Brunello')
    
    # Densities
        # multicrispr::plot_intervals(subset(exonspacers, seqnames==1))
        exonspacers$set <- ifelse(exonspacers$in_brunello, 'brunello', 'multicrispr')
        p_brunello_on <-
            ggplot2::ggplot() + 
            ggplot2::geom_density( data    = data.table::as.data.table(exonspacers), 
                                mapping = ggplot2::aes(x = Doench2016), 
                                fill   = "#00BFC4", alpha = 0.5) + 
            ggplot2::geom_density( data    = data.table::as.data.table(subset(exonspacers, in_brunello==TRUE)), 
                                mapping = ggplot2::aes(x = Doench2016), 
                                fill   = "#F8766D", alpha = 0.5) + 
            ggplot2::theme_bw() + 
            ggplot2::guides(color = FALSE)
        
    # Genome (mis)matches
        color_values <- c(multicrispr = "#00BFC4", brunello = "#F8766D", unique = "#00BA38")
        #autonomics.plot::make_gg_colors(c('brunello', 'multicrispr'), show=TRUE)
        multioff <- table(exonspacers$G0)
        names(multioff) %<>% as.numeric() %>% magrittr::subtract(1) %>% as.character()
        total <- sum(multioff)
        multioff[1] <- 100 * (sum(multioff[1:length(multioff)])) / total
        multioff[2] <- 100 * (sum(multioff[2:length(multioff)])) / total
        multioff[3] <- 100 * (sum(multioff[3:length(multioff)])) / total
        multioff[4] <- 100 * (sum(multioff[4:length(multioff)])) / total
        multioff %<>% extract(1:4)
        
        broff <- table(subset(exonspacers, in_brunello==TRUE)$G0)
        names(broff) %<>% as.numeric() %>% magrittr::subtract(1) %>% as.character()
        total <- sum(broff)
        broff[1] <- 100 * (sum(broff[1:length(broff)])) / total
        broff[2] <- 100 * (sum(broff[2:length(broff)])) / total
        broff[3] <- 100 * (sum(broff[3:length(broff)])) / total
        broff[4] <- 100 * (sum(broff[4:length(broff)])) / total
        broff %<>% extract(1:4)

        dt <- data.table::data.table(
                x   = c('>=1', '>=2', '>=3', '>=1', '>=2', '>=3'), 
                y   = c(multioff[-1], broff[-1]),
                set = c(rep('multicrispr', each=3), rep('brunello', each=3)))
        p_brunello_off<-ggplot2::ggplot(dt) + 
                        ggplot2::geom_col(ggplot2::aes(x=x, y=y, fill =set), position = 'dodge') + 
                        ggplot2::xlab('Offtarget matches') + 
                        ggplot2::ylab('Percentage of spacers') + 
                        ggplot2::theme_bw() + 
                        ggplot2::guides(fill=FALSE)
        
require(multipanelfigure)

fig <- multi_panel_figure(columns = 10, rows = 3, width = 500, height = 250, row_spacing = 10, column_spacing = 10, panel_label_type = 'none')
fig <- fill_panel(fig, p_pe_venn,                row = 1,   column = 1:5 )
fig <- fill_panel(fig, p_pe_on,          row = 2,   column = 1:5 )
fig <- fill_panel(fig, p_parallel_venn,          row = 1,   column = 6:7   )
fig <- fill_panel(fig, p_parallel_on,    row = 2,   column = 6:7 )
fig <- fill_panel(fig, p_brunello_venn,          row = 1,   column = 8:10 )
fig <- fill_panel(fig, p_brunello_on,            row = 2,   column = 8:10 )
fig <- fill_panel(fig, p_brunello_off,           row = 3,   column = 8:10)
fig

require(Cairo)
file_local <- '../multicrisprout/validation.pdf'
file_agnerds <- 'Z:\\abhagwa\\multicrisprmanuscript\\figures\\validation.pdf'
cairo_pdf(file_local, width = figure_width(fig, unit_to = 'inch'), height = figure_height(fig, unit_to = 'inch'))
#pdf(filename, width = figure_width(fig, unit_to = 'inch'), height = figure_height(fig, unit_to = 'inch'))
print(fig)
dev.off()

