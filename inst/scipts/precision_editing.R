# Load 
#======
require(magrittr)
require(multicrispr)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38.masked::BSgenome.Hsapiens.UCSC.hg38.masked
bsinfo <- BSgenome::seqinfo(bsgenome)

# HBB - sickle cell anemia
#=========================

# Sickle cell SNP is at chr11:5227002
# Explore environment around this locus
(fwd <- BSgenome::getSeq(bsgenome, names = 'chr11', start = 5227002-50, end = 5227002+50, strand = '+'))
(rev <- Biostrings::complement(fwd))

# Anzalone et al. (2019): two prime editing sites allowed fixing
# Retrieve sequences (pegspacer, pegext, and nickspacer) to reconstruct precision engineering topology
    # excl needs to be a regex - https://support.bioconductor.org/p/126593/#126606
    excl <- BSgenome::seqnames(bsgenome) %>% setdiff('chr11') %>% paste0('^', ., '$')
    
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
#                                                                                                                                                 ^        3'  nicking spacer 5' 
#                                                                             GGCAGACTTCTC                                                        ^        ACTGAAAATACGGGTCGGGAC        
#                                                                        GTAAC            CaC                                                     ^                             
#                3.5                5'   CTTCATCCACGTTCACCTTGCCCCACAGGGCA                    AGGAGTCAGATGCACCATGGTGTCTGTTTGAGGTTGCTAGTGAACACAG... ^ ...GCCCTGACTTTTATGCCCAGCCCTG   3'
#                                   3'   GAAGTAGGTGCAAGTGGAACGGGGTGTCCCGT                    TCCTCAGTCTACGTGGTACCACAGACAAACTCCAACGATCACTTGTGTC... ^ ...CGGGACTGAAAATACGGGTCGGGAC   5'
#                                                                        CATTGCCGTCTGAAGAGGtG                                                     ^                           
#                                                                                                                                                 ^    5227089           5227111
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
    
    (cas9s <- multicrispr::find_cas9s(gr))

# PRNP snp: Kuru resistance variant (G -> T)
    gr  <-  GenomicRanges::GRanges(
                'chr20:4699500', strand = '+', seqinfo = bsinfo)
    gr  %<>% multicrispr::add_inverse_strand()
    (gr %<>% multicrispr::add_seq(bsgenome))
    (extended <- multicrispr::extend(gr, bsgenome = bsgenome))
    Biostrings::complement(Biostrings::DNAStringSet(extended$seq[2]))
    (cas9s <- multicrispr::find_cas9s(extended))

    # Precision editing
    #
    #    ------------^---===....*....
    # 5' TGGTGGCTGGGG TCAAGGAGGTGGCACCCACAGTCAGTGGAACAA 3'
    # 3' ACCACCGACCCC AGTTCCTCCACCGTGGGTGTCAGTCACCTTGTT 5'
    #

# Tay Sachs    
    gr <- GenomicRanges::GRanges('chr15:72346580-72346583', strand = '-', 
                                 seqinfo = BSgenome::seqinfo(bsgenome))
    
    BSgenome::getSeq(bsgenome, gr)
    gr %<>% multicrispr:::add_inverse_strand()
    extended <- multicrispr::extend(gr)
    extended %<>% multicrispr::add_seq(bsgenome)
    cas9s <- multicrispr::find_cas9s(extended)
    
    # Precision editing
    # 
    #                       ===*================---
    # 5' AATGTGAGACAGCTTAAAATAAAATTAACTATAAGAAACTGGTAA
    # 3' TTACCAGTTTCTTATAGTTAATTTTATTTTAAGCTGTCTCACATT
