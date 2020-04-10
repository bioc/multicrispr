require(magrittr)
require(data.table)
require(multicrispr)

bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38 
x <- multicrispr::char_to_granges(
        c(  PRNP  = 'chr20:4699600:+',  # snp
            HBB  = 'chr11:5227002:-',   # snp
            HEXA = 'chr15:72346579:-'), # ins
        bsgenome = bsgenome)
spacers <- multicrispr::find_pe_spacers(x, bsgenome)
spacers['HBB_1']
    # CATGGTGCATCTGACTCCTG   # HBB_1 C->T, extra nucleotide
    # GCATGGTGCACCTGACTCCTG  # Anzalone
spacers['HEXA_3']
    #  TACCTGAACCGTATATCCTA   # HEXA_3
    # GTACCTGAACCGTATATCCTA  # Anzalone
spacers['PRNP_4']
    # GCAGTGGTGGGGGGCCTTGG # PRNP_4
    # GCAGTGGTGGGGGGCCTTGG # Anzalone

spacers %<>% sort(ignore.strand=TRUE)
spacers$set <- ifelse(names(spacers) %in% c('HBB_1', 'HEXA_3', 'PRNP_4'), 'anzalone', 'multicrispr')

x <- list(multicrispr = names(subset(spacers, targetname=='HBB')), 
        anzalone = names(subset(spacers, targetname=='HBB' & set=='anzalone')))
plot_venn <- function(x, title=''){
    names(x) %<>% paste0(' (', vapply(x, length, integer(1)), ')')
    p <- VennDiagram::venn.diagram(
        x, 
        filename = NULL, 
        fill     = c("#00BFC4", "#F8766D", "#00BA38")[seq_along(x)], 
        col      = c("#00BFC4", "#F8766D", "#00BA38")[seq_along(x)], 
        cat.pos  = rep( 0,   length(x)),
        cat.dist = rep(-0.1, length(x)),
        #cat.col  = c(multicrispr = "#00BFC4", anzalone = "#F8766D"), 
        scaled   = FALSE, 
        main     = title, 
        main.pos = c(0.5, 1.06))
    grid::grid.draw(p)
    p
}

do_plot_venn <- function(selectedtarget){
    grid::grobTree(
        plot_venn(
            list(multicrispr = names(subset(spacers, targetname==selectedtarget)), 
                 Anzalone    = names(subset(spacers, targetname==selectedtarget & set=='anzalone'))),  
            title = selectedtarget))
}

p1 <- grid::grobTree(do_plot_venn('HBB'))
p2 <- grid::grobTree(do_plot_venn('HEXA'))
p3 <- grid::grobTree(do_plot_venn('PRNP'))
p_pe_venn <- gridExtra::grid.arrange(p1, p2, p3, layout_matrix = matrix(1:3, nrow=1))
grid::grid.draw(p_pe_venn)

spacers$in_anzalone <- spacers$set=='anzalone'
multicrispr::plot_intervals(spacers, facet_var = c('seqnames','targetname'), 
                            size_var = 'in_anzalone')

spacers %<>% add_efficiency(bsgenome, 'Doench2016')
ggplot2::ggplot(as.data.table(spacers), ggplot2::aes(x=in_anzalone, y = Doench2016, color=set)) + 
ggplot2::geom_point() + ggplot2::facet_wrap(~ targetname) + ggplot2::theme_bw()


spacers %<>% add_genome_counts(bsgenome)
p_pe_on <-  ggplot2::ggplot(as.data.table(spacers), ggplot2::aes(y=as.factor(as.character(G0-1)), x = Doench2016, color=set)) + 
                    ggplot2::geom_point(size = 4, shape = 15, alpha=0.6) + 
                    ggplot2::facet_wrap(~ targetname) + 
                    ggplot2::theme_bw() + 
                    ggplot2::ylab('Offtarget matches') + 
                    ggplot2::guides(color = FALSE)

spacers$unique <- spacers$G0==1
reticulate::use_condaenv('azienv')
spacers %<>% add_efficiency(bsgenome, method = 'Doench2016')
p <- multicrispr::plot_intervals(
        spacers, facet_var = c('seqnames', 'targetname'), 
        size_var = 'anzalone')


p + ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.5))

p <- multicrispr::plot_intervals(
        spacers, facet_var = c('seqnames', 'targetname'), 
        alpha_var = 'unique', 
        size_var = 'anzalone')
p + ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4))

# DNMT1
    ensdb <- multicrispr::EnsDb.Hsapiens.v98()
    ensdb %<>% filter(UniprotMappingTypeFilter('DIRECT', condition = "=="))
    ensdb %<>% filter(UniprotDbFilter('SWISSPROT', condition = "=="))
    ensembldb::columns(ensdb)
    mycolumns <- c(
        'UNIPROTID', 'PROTEINID',  'TXID',
        'SEQCOORDSYSTEM', 'SEQNAME', 'SEQSTRAND',
        'GENESEQSTART', 'GENESEQEND', 
        'TXSEQSTART', 'TXSEQEND',     'TXCDSSEQSTART', 'TXCDSSEQEND', 
        'PROTEINSEQUENCE')
    dnmt1 <- ensembldb::select(ensdb, 'DNMT1', mycolumns, 'SYMBOL')

x <- char_to_granges('chr19:10133346-10195135:-', bsgenome)
x %<>% multicrispr::add_seq(bsgenome)
x$seq %>% Biostrings::DNAString() %>% Biostrings::reverseComplement() %>% stringi::stri_detect_fixed('GCGGGCTGGAGCTGTTCGCGC')

Biostrings::matchPattern(Biostrings::DNAString('GCGGGCTGGAGCTGTTCGCGC'), bsgenome$chr19, max.mismatch = 2)


x %<>% multicrispr::extend(-200, 200)
x %<>% multicrispr::add_seq(bsgenome)
x$seq %>% stringi::stri_locate_all_fixed('GTACCTGAACCGTATATCCTA')
char_to_granges( GenomicRanges::end(x)+1-187,
                        GenomicRanges::end(x)+1-207, 
                        seqinfo = BSgenome::seqinfo(bsgenome))
BSgenome::getSeq(bsgenome, 'chr15', 72346577, 72346597, strand='-')

AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, 'HEXA', keytype = 'SYMBOL', columns = 'ENTREZID') # 3073
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
AnnotationDbi::keytypes(txbd)
AnnotationDbi::select(txdb, keys = '3073', keytype = 'GENEID', columns = AnnotationDbi::columns(txdb))


GTACCTGAACCGTATATCCTA
ATCCTTCCAGTCAGGGCCAT

BSgenome::getSeq(bsgenome, 'chr15', 72346580, 72346583, strand='+')
BSgenome::getSeq(bsgenome, 'chr15', 72346580, 72346583, strand='-')
