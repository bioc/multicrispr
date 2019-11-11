# Align cas9s to targets and genome
# @param targets   \code{\link[GenomicRanges]{GRanges-class}}
# @param bsgenome  \code{\link[BSgenome]{BSgenome-class}}
# @param outdir string
# @examples
# 
# # Index genome - takes few hours
# bsgenome  <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
# outdir <- file.path(tempdir(), '.autonomics')
# #index_genome(bsgenome, outdir = outdir)
# 
# # Index targets
# bedfile  <- system.file('extdata/SRF.bed', package = 'multicrispr')
# targets0 <- bed_to_granges(bedfile, genome = 'mm10')
# targets  <- double_flank(targets0, leftstart = -200, leftend=-1, rightstart=1, rightend=200)
# targets <- add_seq(targets, bsgenome)
# index_targets(targets, outdir = outdir)
# 
# # Align cas9s to targets and genome
# require(multicrispr)
# cas9s <- find_cas9s(targets)
# outdir <- tempdir()
# name(bsgenome)
# align_cas9s <- function(cas9s, outdir){
#     
#     # Write cas9seqs to fasta
#     uniquecas9s <- Biostrings::DNAStringSet(unique(cas9s$seq))
#     names(uniquecas9s) <- sprintf('cas%d', 1:length(uniquecas9s))
#     
#     cas9fa <- file.path(outdir, 'uniquecas9s.fa')
#     Biostrings::writeXStringSet(uniquecas9s, cas9fa)
#     
#     # Align against targets
#     bowtie_results <- Rbowtie::bowtie(
#                             sequences = cas9fa, 
#                             index     = paste0(outdir, '/targets/targets'), 
#                             f         = TRUE, # fasta input
#                             m         = 10000,
#                             a         = TRUE,
#                             v         = 2, 
#                             norc      = TRUE,
#                             outfile   = paste0(outdir, '/bowtieout.txt'), 
#                             force     = TRUE)
#     bowtie_dt <- data.table::fread(
#         paste0(outdir, '/bowtieout.txt'), 
#         col.names = c('cas9name', 'str', 'target', 'pos', 'cas9seq', 'qual', 
#                       'matches', 'mismatches'))
#     bowtie_dt[ , mismatch := stringi::stri_count_fixed(mismatches, '>')]
#     bowtie_summaries <- bowtie_dt[  , 
#                                     list(   targets0 = sum(mismatch==0), 
#                                             targets1 = sum(mismatch==1), 
#                                             targets2 = sum(mismatch==2)), 
#                                     by = 'cas9seq' ]
#     
#     vcountres <- count_target_matches(uniquecas9s, unique(targetseqs), mismatch = 0)
#     
#     
#     
#     bowtie_results <- Rbowtie::bowtie(
#         sequences = cas9fa, 
#         index     = paste0(outdir, '/targets'), 
#         prefix    = 'targets') 
# }



# @rdname align_cas9s
# index_targets <- function(targets, outdir = '~/.multicrispr'){
#     
#     # Create Names
#     targetfa <- file.path(outdir, 'targets.fa')
#     targetdir  <- targetfa %>% substr(1, nchar(.)-3)
#     
#     # Write to fasta
#     targetseqs <- Biostrings::DNAStringSet(targets$seq)
#     names(targetseqs) <- sprintf('%s:%s-%s(%s)', 
#                                 as.character(GenomicRanges::seqnames(targets)),
#                                 GenomicRanges::start(targets),
#                                 GenomicRanges::end(targets),
#                                 GenomicRanges::strand(targets))
#     Biostrings::writeXStringSet(targetseqs, targetfa)
#     
#     # Index targets
#     Rbowtie::bowtie_build(  targetfa, 
#                             targetdir,  
#                             prefix = basename(targetdir), 
#                             force = TRUE)
#     
# }



# @rdname align_cas9s
# index_genome <- function(bsgenome, outdir = '~/.multicrispr'){
#     
#     # Create Names
#     dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
#     genomefa <- paste0(outdir, '/', bsgenome@pkgname, '.fa')
#     genomedir  <- genomefa %>% substr(1, nchar(.)-3)
#     
#     # Write to fasta
#     BSgenome::writeBSgenomeToFasta(bsgenome, genomefa)
# 
#     # Index genome
#     subdirs <- list.dirs(genomedir, full.names = FALSE, recursive = FALSE)
#     Rbowtie::bowtie_build(  genomefa, 
#                             genomedir, 
#                             prefix = basename(genomedir))
# }


