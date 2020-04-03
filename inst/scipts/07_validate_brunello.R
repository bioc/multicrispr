# Load functions
#===============
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


read_brunello_dt <- function(){
    
    # Read file
    local_brunello_file = '../multicrisprout/brunello/brunello.txt'
    brunello_url <- paste0('https://www.addgene.org/static/cms/filer_public/', 
                           '8b/4c/8b4c89d9-eac1-44b2-bb2f-8fea95672705/', 
                           'broadgpp-brunello-library-contents.txt')
    if (!file.exists(local_brunello_file)){
        download.file(brunello_url, brunello_file)}
    brunello <- data.table::fread(local_brunello_file)
    cmessage('\t%d gRNAs described in Brunello txt file', nrow(brunello))
    
    # Rm positive controls
    n0 <- nrow(brunello)
    brunello %<>% extract(!is.na(`Target Gene ID`))
    ntranscripts <- length(unique(brunello$`Target Transcript`))
    ngenes       <- length(unique(brunello$`Target Gene Symbol`))
    n1 <- nrow(brunello)
    cmessage(paste0('\t%d gRNAs targeting %d transcripts ', 
                    '(%d genes), after removing positive controls'),
             n1, ntranscripts, ngenes)
    
    # Add refseq mrna values
    brunello[, refseq_mrna := `Target Transcript` %>% substr(1, nchar(.)-2) ]
    
    # Return
    brunello
}

load_ensmart <- function(){
    ensembl_dataset <- 'hsapiens_gene_ensembl'
    ensembl_version <- 99
    cmessage('\t\tUse biomaRt: %s version %s', ensembl_dataset, ensembl_version)
    ensmart <- biomaRt::useEnsembl(
                'ensembl', dataset=ensembl_dataset, version=ensembl_version)
    ensmart
}

add_brunello_seqnames <- function(brunello, ensmart){
    
    chromdt <-  suppressMessages(
        biomaRt::getBM( attributes = c('refseq_mrna', 'chromosome_name'), 
                        filters     = 'refseq_mrna', 
                        values      = unique(brunello$refseq_mrna), 
                        mart        = ensmart)) %>% 
        data.table::as.data.table()
    canonical_chromosomes <- GenomeInfoDb::standardChromosomes(
        BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38)
    chromdt %<>% extract(chromosome_name %in% canonical_chromosomes)
    chromdt[, .SD[.N>1], by='refseq_mrna']
    chromdt %>% data.table::setnames('chromosome_name', 'seqnames')
    
    brunello %<>% data.table::merge.data.table(chromdt, by = 'refseq_mrna') 
    cmessage('\t%d gRNAs targeting %d transcripts, after adding canonical chromosome mappings', 
            nrow(brunello), brunello[, length(unique(refseq_mrna))])
    brunello
}

add_brunello_strand <- function(brunello, ensmart){
    
    # Get strand mappings
    brunello %>% data.table::setnames('Strand', 'sense')
    strandt <-  suppressMessages(
        biomaRt::getBM( attributes = c('refseq_mrna', 'strand'), 
                        filters     = 'refseq_mrna', 
                        values      = unique(brunello$refseq_mrna), 
                        mart        = ensmart)) %>%
        data.table::as.data.table()
    strandt %<>% extract(, .SD[.N==1], by = 'refseq_mrna')
    strandt[, strand := vapply(strand, multicrispr::csign, character(1))]
    brunello %<>% data.table::merge.data.table(strandt, by = 'refseq_mrna') 
                    # 76 441 -> 75 463 -> 75 244
    cmessage('\t%d gRNAs targeting %d transcripts, after adding unique strand mappings', 
            nrow(brunello), brunello[, length(unique(refseq_mrna))])
    brunello[sense=='antisense', strand := ifelse(strand=='+', '-', '+')]
    brunello$sense <- NULL
    brunello
}

forge_spacer_ranges <- function(brunello){
    
    # Grange
    brunello %>% data.table::setnames(c('Position of Base After Cut (1-based)'), c( 'start'))
    brunello %>% extract(, end := start)
    bsgenome <- BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
    brunello <- dt2gr(brunello, BSgenome::seqinfo(bsgenome)[GenomeInfoDb::standardChromosomes(bsgenome)])
    
    # Extend
    brunello[GenomicRanges::strand(brunello)=='+'] %<>% extend(-17, +2)
    brunello[GenomicRanges::strand(brunello)=='-'] %<>% extend(-16, +3)
    brunello$crisprspacer <- BSgenome::getSeq(bsgenome, brunello, as.character=TRUE)
    brunello %<>% subset(`sgRNA Target Sequence`==crisprspacer) 
        # Drop 12 non-matches (due to sequence changes). 75 244 -> 75 232
    brunello$`sgRNA Target Sequence` <- NULL
    cmessage('\t%d gRNAs targeting %d transcripts, after removing non-matching spacers', 
            length(brunello), length(unique(brunello$refseq_mrna)))

    # Harmonize/drop names
    names(mcols(brunello)) %<>% 
    stri_replace_first_fixed('Target Context Sequence', 'crisprcontext') %>% 
    stri_replace_first_fixed('PAM Sequence',            'crisprpam')     %>%
    stri_replace_first_fixed('Rule Set 2 score',        'Doench2016')
    brunello$`Genomic Sequence` <- brunello$`Target Gene ID` <- NULL
    brunello$`Exon Number` <- NULL
    brunello$crisprname <- uniquify(brunello$refseq_mrna)
    names(brunello) <- brunello$crisprname
    
    # Return
    brunello
}

reconstruct_brunello_spacers <- function(){
    brunello <- read_brunello_dt()
    ensmart  <- load_ensmart()
    brunello %<>% add_brunello_seqnames(ensmart)
    brunello %<>% add_brunello_strand(ensmart)
    brunello <- forge_spacer_ranges(brunello)
    # Are context seqs identical? Yes!
    #all(brunello$crisprcontext == (brunello %>% extend(-4,+6) %>% getSeq(bsgenome, ., as.character=TRUE)))
    brunello
}

validate_brunello_spacers <- function(brunello){
    
    # Find spacers
    message('Find Brunello spacers')
    bsgenome <- BSgenome::getBSgenome(genome = genome(brunello)[1])
    spacers  <- extend(brunello, 0, +3) %>% 
                multicrispr::find_spacers(bsgenome, plot = FALSE)
        # All brunello spacers are found
        # 4 813 additional spacers are found on rev strands
    
    # Venn 
    spcoords <- unname(as.character(granges(spacers)))
    brcoords <- unname(as.character(granges(brunello)))
    autonomics.plot::plot_venn(list(multicrispr=spcoords, brunello=brcoords))
    spacers %<>% extract(spcoords %in% brcoords) # Limit to Brunello spacers
    
    # Doench score
    message('Doench 2016 score')
    names(mcols(spacers)) %<>% stringi::stri_replace_first_fixed('Doench2016', 'BrunelloDoench2016')
    #brunello_spacers[1:10] %>% add_efficiency(bsgenome, 'Doench2016', plot = FALSE)
    spacers %<>% add_efficiency(bsgenome, 'Doench2016', plot = FALSE, chunksize=10000)
    plot(spacers$BrunelloDoench2016, spacers$Doench2016)
    plot(density(spacers$Doench2016))
    
    # Add genome counts    
    GenomeInfoDb::seqlevelsStyle(spacers) <- 'UCSC'
    bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    spacers %<>% add_genome_counts(bsgenome, mismatches=1)
    table(spacers$G0)
    spacers$unique <- spacers$G0==1
    autonomics.plot::plot_venn(list(
        brunello = spacers$crisprname, 
        unique   = subset(spacers, G0==1)$crisprname))
}

find_enclosing_exons <- function(spacers){
    
    # Load canonical exons
    exons <- ensembldb::exons(multicrispr::EnsDb.Hsapiens.v99())
    exons %<>% extract( as.character(seqnames(exons)) %in% 
                        as.character(GenomicRanges::seqnames(spacers)))
    
    # Find overlaps
    exons   %<>% gr2dt() %>% dt2gr(seqinfo(exons  )[seqlevelsInUse(exons  )])
    spacers %<>% gr2dt() %>% dt2gr(seqinfo(spacers)[seqlevelsInUse(spacers)])
    res <- findOverlaps(exons, extend(spacers, 0, +3), ignore.strand = TRUE, minoverlap = 23)
    
    # Exon spacers
    exonspacers <- spacers %>% extract(unique(subjectHits(res)))
    cmessage('\t%d exon spacers', length(exonspacers))
    
    # Spacer exons
    spacerexons <- exons %>% extract(queryHits(res))
    spacerexons$targetname <- names(spacers)[subjectHits(res)]
    spacerexons %<>% gr2dt()
    spacerexons %<>% extract(, .SD[width == min(width)], by = 'targetname')   # 75 229 / 
    spacerexons %<>% extract(, .SD[1], by = 'targetname')
    spacerexons %<>% dt2gr(seqinfo(exons))
    cmessage('\t%d spacer exons', length(spacerexons))
    
    # Return
    return(spacerexons)

}
    
# Get Brunello exons
#===================
    brunello <- reconstruct_brunello_spacers()
    brunello_spacers <- validate_brunello_spacers(brunello)
    brunello_exons   <- find_enclosing_exons(brunello)
    brunello %<>% extract(brunello_exons$targetname)

# Analyze Brunello exons
#=======================
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
    
# Plot validation results
#========================
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

