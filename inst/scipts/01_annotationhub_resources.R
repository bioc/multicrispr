# Annotation Mus musculus
require(AnnotationHub)
ah <- AnnotationHub()
names(ah)[1:10]
ah <- query(ah, 'Mus musculus')
ah <- query(ah, 'Homo sapiens')
table(ah$dataprovider)

# UCSC
ucsc <- ah[ah$dataprovider=='UCSC']                   # UCSC
table(ucsc$rdataclass)                                #-----
    granges <- query(ucsc, 'GRanges')                 # 173 GRanges with UCSC tracks
    table(granges$sourcetype)
    
    txdbs <- query(ucsc, 'TxDb')                      #   5 TxDb G/T/P models
    table(txdbs)                                      #     mm9/10
    
    twobitfiles <- query(ucsc, 'twoBitFile')          #   4 TwoBitFiles seqs
    twobitfiles                                       #     mm7, 8, 9, 10

# NCBI                                                # NCBI
ncbi <- query(ah, 'NCBI')                             #-----
                                                      #   1 org.Mm.eg.db with gene annotations
    
# Ensembl
ensembl <- ah[ah$dataprovider == 'Ensembl']
table(ensembl$rdataclass)                             # Ensembl
                                                      #---------
    ensdbs <- query(ensembl, 'EnsDb')                 #  12  EnsDb: 87->98
    ensdbs$title
    
    twobitfiles <- query(ensembl, 'twoBitFile')       # 815  TwoBitFile seqs
    sort(table(twobitfiles$genome), decreasing=TRUE)  

    twobitfiles <- query(twobitfiles,'GRCm38')        #  75  GRCm38 (others: BALB_cJ_v1, C3H_HeJ_v1, ...)
    data.frame(names(twobitfiles)[order(twobitfiles$title)], twobitfiles$title[order(twobitfiles$title)])
    twobitfiles['AH49772']                            #       Lots of data redundancy due to periodic automated builds
    twobitfiles['AH50117']                            #       https://support.bioconductor.org/p/125690
    table(twobitfiles$title)                          #       cdna, dna(_r|sm), ncrna

    granges <- query(ensembl, 'GRanges')              # 556 GRanges GTFs
    granges[1]
    sort(granges$title)

