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
spacers$anzalone <- names(spacers) %in% c('HBB_1', 'HEXA_3', 'PRNP_4')
multicrispr::plot_intervals(spacers, facet_var = c('seqnames','targetname'), 
                            size_var = 'anzalone')

spacers %<>% add_genome_counts(bsgenome)
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