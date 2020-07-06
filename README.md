<center> <h1> multicrispr </h1> </center>

[![Build Status](https://travis-ci.com/bhagwataditya/multicrispr.svg?branch=master)](https://travis-ci.com/bhagwataditya/multicrispr)

### Installation

    # Install multicrispr
        remotes::install_git('https://gitlab.gwdg.de/loosolab/software/multicrispr.git', 
        			    repos = BiocManager::repositories())
    # Install azimuth
        install.packages('reticulate')
        reticulate::conda_create('azienv', c('python=2.7'))
        reticulate::use_condaenv('azienv')
        reticulate::py_install(c('azimuth', 'scikit-learn==0.17.1'), 'azienv', pip = TRUE)
    # Index mm10 and hg38
        BiocManager::install('BSgenome.Mmusculus.UCSC.mm10')
        BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')
        index_genome(BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
        index_genome(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)