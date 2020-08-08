require(magrittr)
require(multicrispr)
require(data.table)
require(GenomicRanges)
bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38  
gr <- char_to_granges(c(PRNP = 'chr20:4699600:+',             # snp
                        HBB  = 'chr11:5227002:-',             # snp
                        HEXA = 'chr15:72346580-72346583:-',   # del
                        CFTR = 'chr7:117559593-117559595:+'), # ins
                       bsgenome)
spacers <- find_primespacers(gr, bsgenome)
spacers$targetlocus <- c(rep('chr7 : 117 559 593', 2), 
                        rep('chr11 : 5 227 002',  2),
                        rep('chr15 : 72 346 580',  2),
                        rep('chr20 : 4 699 600',  4))

plot_intervals(spacers, alpha_var = NULL, facet_var = c('targetname', 'targetlocus'), size_var = NULL) + 
ggplot2::guides(color=FALSE, linetype = FALSE, size=FALSE, alpha=FALSE)
ggplot2::ggsave('../graphs/pe_spacers.pdf', width =6, height = 4, device = grDevices::cairo_pdf)

# Matching spacers only (without pam)
(bowtie_results <- multicrispr:::bowtie_count(spacers$crisprspacer, index_genome(bsgenome), norc=FALSE, mismatches = 3))
#(pdict_results  <- multicrispr:::pdict_count( spacers$crisprspacer, bsgenome,               norc=FALSE, mismatches = 3))
#saveRDS(pdict_results, '../offtarget_comparison/pdict_results.rds')
(pdict_results <- readRDS('../offtarget_comparison/pdict_results.rds'))


# spacer+pam (expanded)
spacers %<>% extract(.$targetname=='PRNP')
mcols(spacers) %<>% extract(, 1:10)
(bowtie_results2 <- multicrispr:::count_offtargets(spacers[2], bsgenome, mismatches=1, offtargetmethod = 'bowtie'))
(pdict_results2  <- multicrispr:::count_offtargets(spacers[2], bsgenome, mismatches=1, offtargetmethod = 'pdict' ))


bowtie_results[c("GCAGCTGGGGCAGTGGTGGG"), on = 'readseq']
pdict_results[c("GCAGCTGGGGCAGTGGTGGG"),  on = 'readseq']
mcols(bowtie_results2)[, c('crisprspacer', 'off0', 'off1')]
mcols(pdict_results2)[, c('crisprspacer', 'off0', 'off1')]


    




spacers2 <- spacers %>% vadd_genome_matches(bsgenome, mismatches=1, verbose=TRUE)




merge(spacerdt, gr2dt(vres), by = 'index')





(pdict_results2  <- multicrispr:::count_spacer_matches(spacers,              bsgenome,  norc=FALSE, mismatches = 3, offtargetmethod = 'vcountpdict'))

offtarget_comparison <- data.table( targetname   = spacers$targetname, 
            crisprname   = spacers$crisprname,
            crisprrange  = as.character(granges(spacers)),
            spacer       = spacers$crisprspacer,
            pam          = spacers$crisprpam, 
            V0           = pdict_results[ spacers$crisprspacer, on = 'readseq']$MM0,
            B0           = bowtie_results[spacers$crisprspacer, on = 'readseq']$MM0, 
            V1           = bowtie_results[spacers$crisprspacer, on = 'readseq']$MM1,
            B1           = pdict_results[ spacers$crisprspacer, on = 'readseq']$MM1,
            V2           = pdict_results[ spacers$crisprspacer, on = 'readseq']$MM2,
            B2           = bowtie_results[spacers$crisprspacer, on = 'readseq']$MM2,
            V3           = pdict_results[ spacers$crisprspacer, on = 'readseq']$MM3,
            B3           = bowtie_results[spacers$crisprspacer, on = 'readseq']$MM3
)
fwrite(offtarget_comparison, '../offtarget_comparison/offtarget_comparison.txt', sep = '\t')
