# multicrispr: gRNA design

[![](https://bioconductor.org/shields/build/devel/bioc/multicrispr.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/multicrispr/)
[![](https://bioconductor.org/shields/availability/3.12/multicrispr.svg)](https://bioconductor.org/packages/devel/bioc/html/multicrispr.html#archives) 
[![](https://bioconductor.org/shields/years-in-bioc/multicrispr.svg)](https://bioconductor.org/packages/devel/bioc/html/multicrispr.html#since)
[![](https://img.shields.io/badge/doi-10.26508/lsa.202000757-blue.svg)](https://doi.org/10.26508/lsa.202000757)


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


### Overview
   
![](https://gitlab.gwdg.de/loosolab/software/multicrispr/-/wikis/uploads/cdf31586bcf776a7a40acaaaf5172e10/overview.png)

