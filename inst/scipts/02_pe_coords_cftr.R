# Load functions
    require(magrittr)
    require(multicrispr)
    filter <- ensembldb::filter
    UniprotMappingTypeFilter <- ensembldb::UniprotMappingTypeFilter
    UniprotDbFilter          <- ensembldb::UniprotDbFilter

    
# Map identifiers: CFTR
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
    cftr <- ensembldb::select(ensdb, 'CFTR', mycolumns, 'SYMBOL')
    cftr # Two transcripts
         # Two distinct proteins
         # ENSP00000003084 gold labeled
    ensp <- 'ENSP00000003084'
    cftr %<>% extract(.$PROTEINID == ensp, )
    cftr$PROTEINSEQUENCE %>% substr(507, 509) # DELTA F 508
    gr   <- ensembldb::proteinToGenome(
                IRanges::IRanges(start=508, end=508) %>% set_names(ensp), 
                ensdb) %>% 
            extract2(ensp) %>% 
            (function(y){seqlevelsStyle(y) <- 'UCSC'; y})
    bs <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    BSgenome::getSeq(
        bs, seqnames(gr), start = start(gr), end = end(gr), strand = '+')
