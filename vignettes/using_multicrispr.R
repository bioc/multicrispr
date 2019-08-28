## ----setup, include = FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "    #", 
  cache = TRUE
)

## ----read------------------------------------------------------------------
library(magrittr)
library(multicrispr)
bedfile      <- system.file('extdata/SRF.bed', package = 'multicrispr')
bsgenome     <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
targetranges <- read_bed(bedfile, bsgenome) %>% flank_fourways()

