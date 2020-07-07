#------------------------------------------------------------------------
# Protein coordinates -> Genomic coordinates
#              HBB: E6V : GGC>GTC : chr20:4699600
#------------------------------------------------------------------------

# Load
    require(magrittr)
    require(multicrispr)
    filter <- ensembldb::filter
    UniprotMappingTypeFilter <- ensembldb::UniprotMappingTypeFilter
    UniprotDbFilter          <- ensembldb::UniprotDbFilter

    bs <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    ensdb <- multicrispr::EnsDb.Hsapiens.v98()
    ensdb %<>% filter(UniprotMappingTypeFilter('DIRECT', condition = "=="))
    ensdb %<>% filter(UniprotDbFilter(      'SWISSPROT', condition = "=="))

# Map identifiers: HBB -> ENST00000335295
    ensembldb::columns(ensdb)
    mycolumns <- c(
        'UNIPROTID', 'PROTEINID',  'TXID',
        'SEQCOORDSYSTEM', 'SEQNAME', 'SEQSTRAND',
        'GENESEQSTART', 'GENESEQEND', 
        'TXSEQSTART', 'TXSEQEND',     'TXCDSSEQSTART', 'TXCDSSEQEND', 
        'PROTEINSEQUENCE')
    hbb <- ensembldb::select(ensdb, 'HBB', mycolumns, 'SYMBOL')
            # Two transcript isoforms: ENST00000335295  &  ENST00000647020
            # Translate into single proteform, differences are in 3'UTR
            # Both are gold labeled in ensembl browser, T*95 has a RefSeq match
            # Let's use ENST00000335295, since it has a RefSeq match
            # (For genomic mapping it actually doesnt matter)
    hbb %<>% extract(.$TXID=='ENST00000335295', )
    
# Map coordinates: E6V -> GAG (codes for E) -> chr11:5227002
    ensp <- 'ENSP00000333994'
    substr(hbb$PROTEINSEQUENCE[1], 1, 7)
        # It's E7, not E6
        # Reason: initiator Methionine is cleaved off!
        # Mature protein doesn't contain it.
        # However, Uniprot Sequence does contain it!
    gr   <- ensembldb::proteinToGenome(
                IRanges::IRanges(start=7, end=7) %>% set_names(ensp), 
                ensdb) %>% 
            extract2(ensp) %>% 
            (function(y){seqlevelsStyle(y) <- 'UCSC'; y})
    BSgenome::getSeq(
        bs, seqnames(gr), start = start(gr), end = end(gr), strand = '-')
    gr
    
# Explore environment around this locus
bsgenome <- BSgenome.Hsapiens.UCSC.hg38.masked::BSgenome.Hsapiens.UCSC.hg38.masked
bsinfo <- seqinfo(bsgenome)
(fwd <- BSgenome::getSeq(bsgenome, names = 'chr11', start = 5227002-50, end = 5227002+50, strand = '+'))
(rev <- Biostrings::complement(fwd))

# Anzalone et al. (2019): two prime editing sites allowed fixing
# Retrieve sequences (pegspacer, pegext, and nickspacer) to reconstruct precision engineering topology
    # excl needs to be a regex - https://support.bioconductor.org/p/126593/#126606
    excl <- seqnames(bsgenome) %>% setdiff('chr11') %>% paste0('^', ., '$')
    
    # 3.5
    (pegspacer3.5   <- BSgenome::vmatchPattern('GTAACGGCAGACTTCTCCAC',        bsgenome, min.mismatch = 0, max.mismatch = 2, exclude = excl))
    (pegext3.5      <- BSgenome::vmatchPattern('ACCTGACTCCTGAGGAGAAGTCTGCC',  bsgenome, min.mismatch = 0, max.mismatch = 2, exclude = excl))
    (nickspacer3.5  <- BSgenome::vmatchPattern('GGGCTGGGCATAAAAGTCA',         bsgenome, min.mismatch = 0, max.mismatch = 2, exclude = excl))
    
    # 3.7
    (pegspacer3.7   <- BSgenome::vmatchPattern('GCATGGTGCACCTGACTCCTG',       bsgenome, min.mismatch = 0, max.mismatch = 2, exclude = excl))
    (pegext3.7      <- BSgenome::vmatchPattern('AGACTTCTCCTCAGGAGTCAGGTGCAC', bsgenome, min.mismatch = 0, max.mismatch = 2, exclude = excl))
    (nickspacer3.7  <- BSgenome::vmatchPattern('GCCTTGATACCAACCTGCCCA',       bsgenome, min.mismatch = 0, max.mismatch = 2, exclude = excl))

# Reconstruct prime editing topology
# Note: GG is on the non-targeted strand
#       H840A nickase used, which cleaves non-target strand.
#                                                                                primer    f       RT   
#                                                                             bindingsite  i    template
#                                                                                          x          C
#                                                                          3' CCGTCTGAAGAGGaGTCCTCAGTCTA 5'<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
#                                                                             ::::::::::::                                                        ^        3'  nicking spacer 5' 
#                                                                             GGCAGACTTCTC                                                        ^        ACTGAAAATACGGGTCGGGAC        
#                                                                        GTAAC            CaC                                                     ^                             
#                3.5                5'   CTTCATCCACGTTCACCTTGCCCCACAGGGCA                    AGGAGTCAGATGCACCATGGTGTCTGTTTGAGGTTGCTAGTGAACACAG... ^ ...GCCCTGACTTTTATGCCCAGCCCTG   3'
#                                   3'   GAAGTAGGTGCAAGTGGAACGGGGTGTCCCGT                    TCCTCAGTCTACGTGGTACCACAGACAAACTCCAACGATCACTTGTGTC... ^ ...CGGGACTGAAAATACGGGTCGGGAC   5'
#                                                                        CATTGCCGTCTGAAGAGGtG                                                     ^                           
#                                                                        ::::::::::::::::::::                                                     ^    5227089           5227111
#                                                                    5'  GTAACGGCAGACTTCTCCaC >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>                                    
#                                                                                                               
#                                                                           pegRNAspacer         
#            
#                
#                                                                                               pegRNAspacer
#                                                             pegRNA                                          C         G
#                3.7                 V<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 3' GTCCTCAGTCTACGTGGTACC  5'
#                                    V                                                                                 
#                                    V                                                      CAGGAGTCAGATGCACCATGG
#   5'  ACCTTGATACCAACCTGCCCACCCG... V ..CTTCATCCACGTTCACCTTGCCCCACAGGGCAGTAACGGCAGACTTCTCCa                     TGTCTGTTTGAGGTTGCTAGTGAACACAG     3'
#   3'  TGGAACTATGGTTGGACGGGTGGGC... V ..GAAGTAGGTGCAAGTGGAACGGGGTGTCCCGTCATTGCCGTCTGAAGAGGt                     ACAGACAAACTCCAACGATCACTTGTGTC     5'
#                                    V                                                      GTC             GTACC
#   5'  ACCTTGATACCAACCTGCCCACCCG 3' V                                                         CTCAGTCTACGTG
#        G                           V                                                                          
#              nickingspacer         >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> 5' AGACTTCTCCaCAGGAGTCAGGTGCAC  3'
#                                                                                          f                
#                                                                                   RT     i     primer       
#                                                                                template  x   bindingsite
#
