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
spacers
bowtie_results <- multicrispr:::bowtie_count(spacers$crisprspacer, index_genome(bsgenome), norc=FALSE, mismatches = 3)
pdict_results  <- multicrispr:::pdict_count( spacers$crisprspacer, bsgenome,               norc=FALSE, mismatches = 3)

saveRDS(pdict_results, '../offtarget_comparison/pdict_results.rds')

bowtie_results0 <- bowtie_results
pdict_results0  <- pdict_results

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