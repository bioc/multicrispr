## ---- include=TRUE, echo=FALSE, message=FALSE---------------------------------
knitr::opts_chunk$set(echo = TRUE, collapse = TRUE, cache = TRUE)
#str(knitr::opts_chunk$get())

## ---- fig.small = TRUE, echo = FALSE, fig.cap = 'gRNAs for Crispr/Cas9 (A) and Prime Editing'----
knitr::include_graphics("../inst/extdata/pe.png")

## ----overview, fig.wide = TRUE, out.width = "80%", echo = FALSE---------------
knitr::include_graphics("../inst/extdata/readme_portrait.png")

## ---- eval = FALSE------------------------------------------------------------
#  url <- 'https://gitlab.gwdg.de/loosolab/software/multicrispr.git'
#  remotes::install_git(url, repos = BiocManager::repositories())

## ---- eval = FALSE------------------------------------------------------------
#  # Install - run R(Studio) with admin privileges for this to work!
#    reticulate::conda_create('azienv', 'python=2.7')
#    reticulate::conda_install('azienv', 'azimuth', pip = TRUE)
#    reticulate::conda_install('azienv', 'scikit-learn==0.17.1', pip = TRUE)

## -----------------------------------------------------------------------------
# Activate
  reticulate::use_condaenv('azienv')

## ---- eval = FALSE------------------------------------------------------------
#  index_genome(BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
#  index_genome(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  )

## ---- echo = FALSE, results = FALSE, message=FALSE----------------------------
  # Not required
  # Done to load dependencies silently - keeping focus
  require(GenomicRanges)
  require(Biostrings)
  require(dplyr)
  require(dbplyr)
  require(htmltools)
  require(htmlwidgets)

## ---- fig.small = TRUE, out.width = "70%", echo = FALSE-----------------------
knitr::include_graphics("../inst/extdata/01_define_targets.png")

## -----------------------------------------------------------------------------
require(magrittr)
require(multicrispr)
bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
tfbs0 <- bed_to_granges(bedfile, genome = 'mm10')

## -----------------------------------------------------------------------------
require(multicrispr)
entrezfile <- system.file('extdata/SRF.entrez', package = 'multicrispr')
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
sometfbs0 <- genefile_to_granges(entrezfile, txdb, complement = TRUE)

## ---- fig.width=3.5, fig.height=1.5-------------------------------------------
# char_to_granges: Anzalone et al. (2019) prime editing targets
bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
pe0 <- char_to_granges(c(HBB  = 'chr11:5227002:-'), bsgenome)
plot_intervals(pe0)

## ---- fig.small = TRUE, out.width = "70%", echo = FALSE-----------------------
knitr::include_graphics("../inst/extdata/02_transform.png")

## ---- fig.width=4, fig.height=2, out.width="65%"------------------------------
# Extend
targets <- pe0
invisible(extend(targets, -22, 22, plot = TRUE))
# Up flank
invisible(up_flank(  targets, -200, -1, plot = TRUE))
# Down flank
invisible(down_flank( targets, 1, 200, plot = TRUE))
# Double flank
invisible(double_flank(targets, -200, -1, +1, +200, plot = TRUE))

## ---- fig.small = TRUE, out.width = "70%", echo = FALSE-----------------------
knitr::include_graphics("../inst/extdata/03_find.png")

## ---- fig.width=3, fig.height=1.5---------------------------------------------
bsgenome   <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
neighbourhood <- extend_for_pe(targets, bsgenome)
spacers    <- find_spacers(neighbourhood, bsgenome, complement = FALSE)

## ---- fig.width=3.3, fig.height=1.7-------------------------------------------
spacers <- find_pe_spacers(targets, bsgenome=bsgenome) 

## ---- fig.width=3.6, fig.height=2.1-------------------------------------------
spacers <- find_pe_spacers(targets, bsgenome=bsgenome, nrt = 48) 

## ---- fig.small = TRUE, out.width = "70%", echo = FALSE-----------------------
knitr::include_graphics("../inst/extdata/04_offtargets.png")

## ---- fig.width=3.6, fig.height=2.1-------------------------------------------
if (has_been_indexed(bsgenome)){
  spacers %<>% add_offtargets(bsgenome, mismatch = 0, plot = TRUE)
}

## ---- fig.small = TRUE, out.width = "70%", echo = FALSE-----------------------
knitr::include_graphics("../inst/extdata/05_efficiency.png")

## ---- fig.width=3.6, fig.height=2.1-------------------------------------------
if (reticulate::py_module_available('azimuth')){
  spacers %<>% add_efficiency(bsgenome, 'Doench2016')
}

## -----------------------------------------------------------------------------
  spacers

