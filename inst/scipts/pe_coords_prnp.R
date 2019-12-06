#------------------------------------------------------------------------
# Protein coordinates -> Genomic coordinates
#              PRNP: G127V : GGC>GTC : chr20:4699600
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
    prnp <- ensembldb::select(ensdb, 'PRNP', mycolumns, 'SYMBOL')
    prnp$PROTEINSEQUENCE[[1]] == prnp$PROTEINSEQUENCE[[2]] # identical!
        # Conclusion: ENST00000379440 & ENST00000430350
        #    - differ in length of 3' UTR and 5'UTR
        #    - have identical codingtranscript and protein
        # 
        # For protein coordinate -> genome coordinate mapping
        #    - both are equivalent 
        #    - let's choose ENST00000379440, since it is gold colored in Ensembl 
        #      genome browser (manual curation and automated annotation match)

    
# Map coordinates: G127V  :  GGC>GTC  :  chr20:4699600
    ensp <- 'ENSP00000368752'
    gr   <- ensembldb::proteinToGenome(
                IRanges::IRanges(start=127, end=127) %>% set_names(ensp), 
                ensdb) %>% 
            extract2(ensp) %>% 
            (function(y){seqlevelsStyle(y) <- 'UCSC'; y})
    bs <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    BSgenome::getSeq(bs, names = seqnames(gr), start = start(gr), end = end(gr))
