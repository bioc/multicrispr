# Shariati paper
#    1. Target oct4 binding site 137-151 nt upstream of Nanog TSS
#    2. Target Utf1 gene    

# Nanog: three isoforms, each with a different TSS
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
AnnotationDbi::keytypes(txdb)
ensdb <- multicrispr::EnsDb.Mmusculus.v98()
ensembldb::keytypes(ensdb)
ensembldb::columns(ensdb)
ensembldb::select(ensdb, 'Nanog', keytype = 'GENENAME', columns = ensembldb::columns(ensdb))
ensembldb::select(ensdb, 'Nanog', keytype = 'GENENAME', columns = c('SEQNAME', 'SEQSTRAND', 'TXID', 'TXSEQSTART', 'TXCDSSEQSTART', 'TXCDSSEQEND', 'TXSEQEND', 'GENESEQSTART', 'GENESEQEND'))

# Find spacers in oct4 site upstream of nanog tss
targets <- multicrispr::char_to_granges(c(oct4_nanog='chr6:122707565:+'), bsgenome) %>%
        multicrispr::up_flank(-151, -137)
targets$targetstart <- start(targets)
targets$targetend <- end(targets)
targets %<>% multicrispr::extend(-22, +22)
bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
pattern <- 'GGTCACCTTACAGCTTCTTTTGCATTACAATGTCCATGGTGG'
BSgenome::getSeq(bsgenome, targets, as.character=TRUE) %>% 
stringi::stri_detect_fixed('GGTCACCTTACAGCTTCTTTTGCATTACAATGTCCATGGTGG')
spacers <- multicrispr::find_spacers(targets, bsgenome)
spacers %<>% multicrispr:::add_genome_counts(bsgenome)
reticulate::use_condaenv('azienv')
spacers %<>% multicrispr::add_efficiency(bsgenome, method = 'Doench2016')

# Find spacers in oct4 site downstream of Utf1
ensembldb::select(ensdb, 'Utf1', keytype = 'GENENAME', columns = c('SEQNAME', 'SEQSTRAND', 'TXID', 'TXSEQSTART', 'TXCDSSEQSTART', 'TXCDSSEQEND', 'TXSEQEND', 'GENESEQSTART', 'GENESEQEND'))
targets <- multicrispr::char_to_granges(c(oct4_utf1='chr7:139943789:+'), bsgenome) %>%
        multicrispr::down_flank(+1825, +1825)
targets$targetstart <- start(targets)
targets$targetend <- end(targets)
targets %<>% multicrispr::extend(-100, +100)
BSgenome::getSeq(bsgenome, targets) %>% 
#Biostrings::reverseComplement() %>% 
as.character() %>% 
stringi::stri_detect_fixed('GTTGTTATGCTAGTGAAGTGC')
