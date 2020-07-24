# General 
#========
    # Install azimuth for scoring
    #----------------------------
    # require(reticulate)
    # conda_create('azienv', c('python=2.7'))
    # use_condaenv('azienv')
    # py_install(c('azimuth', 'scikit-learn==0.17.1'), 'azienv', pip = TRUE)
    # reticulate::use_condaenv('azienv')

    # Index genome
    #-------------
    require(multicrispr)
    index_genome(BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
    index_genome(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38) 
    
    # Load packages
    require(magrittr)
    require(multicrispr)
    reticulate::use_condaenv('azienv')
        
# Parallel Targeting
#===================
    
    # Define targets to block TFBS
    bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
    targets  <- bed_to_granges(bedfile, 'mm10') %>% 
                extend(-22, +22)
            
    # Find specific, efficient spacers
    bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    spacers <-  find_spacers(targets, bsgenome)
    spacers %<>% filter_target_specific(targets, bsgenome)
    spacers %<>% add_efficiency(bsgenome, method = 'Doench2016')
    spacers
    
# Prime Editing
#==============
    
    # Define target for editing
    bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
    targets  <- char_to_granges(c(  PRNP = 'chr20:4699600:+',             # snp
                                    HBB  = 'chr11:5227002:-',              # snp
                                    HEXA = 'chr15:72346580-72346583:-',    # del
                                    CFTR = 'chr7:117559593-117559595:+'),  # ins
                                bsgenome)
    
    # Find specific, efficient spacers
    spacers <-  find_primespacers(targets, bsgenome)
    spacers %<>% filter_prime_specific(bsgenome)
    spacers %<>% add_efficiency(bsgenome, method = 'Doench2016')
    spacers
    
