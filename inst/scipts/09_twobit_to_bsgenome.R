require(magrittr)
# AnnotationHub has 2bit files for 224 organisms
ah <- AnnotationHub::AnnotationHub()
ah %<>% extract(.$rdataclass == 'TwoBitFile')
length(unique(ah$species))
table(ah$species)

# Create a BSgenome for SarsCov2
require(BSgenome)
url <- 'http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/bigZips/wuhCor1.2bit'
download.file(
    url, paste0('../multicrisprout/indexedgenomes/SarsCov2', basename(url)))
forgeBSgenomeDataPkg(
    '../multicrisprout/indexedgenomes/SarsCov2/BSgenome.SarsCov2.UCSC.wuhCor1-seed', 
    destdir = '../multicrisprout/indexedgenomes/SarsCov2/SarsCov2')
bspackage <- '../multicrisprout/indexedgenomes/SarsCov2/SarsCov2/BSgenome.SarsCov2.UCSC.wuhCor1'
devtools::check(bspackage)
devtools::install(bspackage)
require(BSgenome.SarsCov2.UCSC.wuhCor1)
bsgenome <- BSgenome.SarsCov2.UCSC.wuhCor1::BSgenome.SarsCov2.UCSC.wuhCor1
#devtools::build(bspackage)

# Multicrispr
# Suppose we want to study the role of a SarsCov2 protein through a Crispr-based
# technology. One such protein is the protein ORF3, for which we see that the 
# UCSC genome browser currently has no known function
require(multicrispr)
gr <- char_to_granges('NC_045512v2:25393-26220:+', bsgenome)
BSgenome::getSeq(bsgenome, gr)
spacers <- multicrispr::find_spacers(gr, bsgenome, complement = FALSE)
index_genome(bsgenome)
spacers %<>% add_genome_counts(bsgenome)
    # Is this correct?

# Suppose we 
spacers %<>%
nsp5 <- 'agugguuuuagaaaaauggcauucccaucugguaaaguugaggguuguaugguacaaguaacuugugguacaacuacacuuaacggucuuuggcuugaugacguaguuuacuguccaagacaugugaucugcaccucugaagacaugcuuaacccuaauuaugaagauuuacucauucguaagucuaaucauaauuucuugguacaggcugguaauguucaacucaggguuauuggacauucuaugcaaaauuguguacuuaagcuuaagguugauacagccaauccuaagacaccuaaguauaaguuuguucgcauucaaccaggacagacuuuuucaguguuagcuuguuacaaugguucaccaucugguguuuaccaaugugcuaugaggcccaauuucacuauuaaggguucauuccuuaaugguucaugugguaguguugguuuuaacauagauuaugacugugucucuuuuuguuacaugcaccauauggaauuaccaacuggaguucaugcuggcacagacuuagaagguaacuuuuauggaccuuuuguugacaggcaaacagcacaagcagcugguacggacacaacuauuacaguuaauguuuuagcuugguuguacgcugcuguuauaaauggagacaggugguuucucaaucgauuuaccacaacucuuaaugacuuuaaccuuguggcuaugaaguacaauuaugaaccucuaacacaagaccauguugacauacuaggaccucuuucugcucaaacuggaauugccguuuuagauaugugugcuucauuaaaagaauuacugcaaaaugguaugaauggacguaccauauuggguagugcuuuauuagaagaugaauuuacaccuuuugauguuguuagacaaugcucagguguuacuuuccaa'
nsp5 %<>% toupper()
nsp5 %<>% stringi::stri_replace_all_fixed('U', 'T')
stringi::stri_detect_fixed(bsgenome$NC_045512v2, nsp5)
bsgenome$NC_045512v2