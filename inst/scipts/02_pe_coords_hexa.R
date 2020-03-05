#------------------------------------------------------------------------
# Protein coordinates -> Genomic coordinates
#              HEXA: 
#------------------------------------------------------------------------

# Load functions
    require(magrittr)
    require(multicrispr)
    filter <- ensembldb::filter
    UniprotMappingTypeFilter <- ensembldb::UniprotMappingTypeFilter
    UniprotDbFilter          <- ensembldb::UniprotDbFilter

    
# Map identifiers: PRNP -> ENST00000379440
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
    hexa <- ensembldb::select(ensdb, 'HEXA', mycolumns, 'SYMBOL')
        # single ENST / ENSP

# Google helped to locate Tay Sachs locus
    gr <- GenomicRanges::GRanges('chr15:72346580-72346583', strand = '-', 
                                 seqinfo = seqinfo(bsgenome))
    BSgenome::getSeq(bsgenome, gr)
