require(magrittr)
require(data.table)
require(multicrispr)
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

# Find spacers in oct4 site ups <- <- <- me, targets) %>% 
#Biostrings::reverseComplement() %>% 
as.character() %>% 
stringi::stri_detect_fixed('GTTGTTATGCTAGTGAAGTGC')
