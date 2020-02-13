# General 
#========
    # Install azimuth for scoring
    #----------------------------
    # require(reticulate)
    # conda_create('azienv', c('python=2.7'))
    # use_condaenv('azienv')
    # py_install(c('azimuth', 'scikit-learn==0.17.1'), 'azienv', pip = TRUE)
    # reticulate::use_condaenv('azienv')

    # Load packages and genome. Index genome.
    #----------------------------------------
    require(magrittr)
    require(multicrispr)
    require(BSgenome.Mmusculus.UCSC.mm10) 
    require(BSgenome.Hsapiens.UCSC.hg38)  # human genome
    indexedhuman <- index_genome(BSgenome.Mmusculus.UCSC.mm10)
    indexedmouse <- index_genome(BSgenome.Hsapiens.UCSC.hg38) 
    
# Parallel Targeting
#===================
    # Define targets for blocking
    bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
    targets  <- bed_to_granges(bedfile, 'mm10') %>% 
                extend(-22, +22)
            
    # Find specific, efficient spacers
    spacers <-  targets %>% 
                find_spacers(bsgenome) %>% 
                filter_target_specific(targets, bsgenome, indexedgenome) %>% 
                add_efficiency(method = 'Doench2016', condaenv = 'azienv')
            