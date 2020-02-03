require(magrittr)
require(multicrispr)
bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
targets <- bed_to_granges(bedfile, genome='mm10', plot = FALSE)
targets %<>% extend()
bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
spacers <- find_spacers(targets, bsgenome, plot = FALSE)
indexedgenome <- "~/.multicrispr/bowtie/BSgenome.Mmusculus.UCSC.mm10"
    # Change this to path of bowtie indexed genome on your system
spacers %<>% add_genome_counts(indexedgenome, plot = FALSE)

#'         # conda create --name azimuthenv python=2.7
#'         # conda activate azimuthenv
#'         # pip install azimuth
#'         # pip install scikit-learn==0.17.1
spacers %<>% add_efficiency(
                bsgenome, method = 'Doench2016', condaenv = 'azimuthenv')
