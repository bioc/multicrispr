# Load functions
require(AnnotationDbi)
require(BSgenome)
require(multicrispr)
require(magrittr)
require(biomaRt)
require(data.table)

# Load Brunello library
brunello_url <- 'https://www.addgene.org/static/cms/filer_public/8b/4c/8b4c89d9-eac1-44b2-bb2f-8fea95672705/broadgpp-brunello-library-contents.txt'
brunello <- data.table::fread(brunello_url)
brunello[, .(N = length(unique(`Target Transcript`))), by = 'Target Gene ID'][, table(N)] #   1 transcript per gene
brunello[, .(N = length(unique(`Exon Number`))),       by = 'Target Gene ID'][, table(N)] # 1-4 exons      per gene
brunello %<>% extract(!is.na(`Target Gene ID`))

# Keep ensdb genes
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene 
ensdb <- multicrispr::EnsDb.Hsapiens.v98()
autonomics.plot::plot_venn(list(brunello = unique(as.character(brunello$`Target Gene ID`)), 
                                txdb     = unique(keys(txdb, keytype = 'GENEID')), 
                                ensdb    = unique(as.character(keys(ensdb, keytype = 'ENTREZID')))))
# Add seqnames
brunello %>% setkey('Target Gene ID')
brunello %<>% extract(`Target Gene ID` %in% keys(ensdb, keytype = 'ENTREZID'))
chroms <- select(ensdb, as.character(brunello$`Target Gene ID`), keytype = 'ENTREZID', columns = 'SEQNAME')
chroms %<>% data.table::as.data.table()
chroms[, table(SEQNAME)]
chroms %<>% extract( , .SD[nchar(SEQNAME) == min(nchar(SEQNAME))], by = 'ENTREZID')
chroms[, table(SEQNAME)]
chroms[,  .SD[.N>1],by = 'ENTREZID' ]
chroms %<>% extract(!(ENTREZID == 266 & SEQNAME=='Y'))
brunello %<>% merge(chroms, by.x = 'Target Gene ID', by.y = 'ENTREZID', all.x=TRUE, sort = FALSE)

# Overlap with exons
brunello %>% setnames('Position of Base After Cut (1-based)', 'start')
brunello[, end:=start]
brunello %>% setnames('SEQNAME', 'seqnames')
brunello$Strand <- '*'
brunellogr <- brunello %>% as('GRanges')
exons1 <- exons(ensdb)

exons_in_brunello <- unique(queryHits(findOverlaps(exons1, brunellogr)))
brunello_in_exons <- unique(subjectHits(findOverlaps(ensex, brunellogr)))
exons1[exons_in_brunello]
brunellogr[brunello_in_exons]
brunello_not_in_exons <- seq_along(brunellogr) %>% setdiff(brunello_in_exons)
brunellogr %<>% extract(-brunello_not_in_exons)
exons1 %<>% extract(exons_in_brunello)

# Find crispr in brunello exons
bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
seqlevelsStyle(bsgenome) <- 'NCBI'
dt <- gr2dt(exons1)
dt$seqnames %<>% droplevels()
exons1 <- dt2gr(dt, seqinfo(exons1)[levels(droplevels(seqnames(exons1)))])
genome(exons1)[] <- 'hg38'
spacers <- multicrispr::find_spacers(exons1, bsgenome, plot = FALSE)
spacers %<>% add_genome_counts(bsgenome, )

brunello #    76 065
spacers  # 4 374 614
spacers %<>% add_genome_counts(mismatches=1)


seqinfo(exons1) <- seqinfo(bsgenome)
bsgenome

columns(tx)
select(tx, keys='1', columns = c('CDSCHROM', 'CDSSTRAND'), keytype = 'GENEID')
getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 'chr19', 58351502-3, 58351502+16, strand='-', as.character = TRUE) %>%
equals(brunello[1]$`sgRNA Target Sequence`)

# one transcript per gene is hit

# Keep the part of brunello that tallies with current ensdb
brunello %>% setkey('Target Gene ID')
brunello[unique(keys(ensdb, keytype = 'ENTREZID'))]


# Map RefSeq mRNA ids to ensembl transcript ids (which can then be mapped to exon rank)
'NM_130786.3'
ensmart <- useMart('ensembl')
listDatasets(ensmart) %>% subset(dataset %>% stringi::stri_detect_fixed('hsapiens'))
ensmart <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
listAttributes(ensmart)
listAttributes(ensmart) %>% extract2('name') %>% extract(stringi::stri_detect_fixed(., 'exon'))
values <- brunello[, `Target Transcript` %>% tstrsplit('[.]') %>% extract2(1)] %>% unique() %>% na.exclude()
bmdt <- getBM(attributes = c('refseq_mrna', 'ensembl_transcript_id'), filters = 'refseq_mrna', values = values, mart = ensmart)
bmdt %<>% as.data.table()
bmdt[, .(N=length(unique(ensembl_transcript_id))), by = 'refseq_mrna'][, table(N)]

brunello$`Target Gene ID`
brunello %<>% extract(!is.na(`Target Gene ID`))

keys1 <- unique(as.character(brunello$`Target Gene ID`))
keys1 %<>% extract(!is.na(keys1))

# Map chr using txdb
txgenes <- select(tx, keys=,keys1, columns = c('CDSCHROM'), keytype = 'GENEID')
txgenes %<>% as.data.table()
txgenes %<>% extract(, .SD[nchar(CDSCHROM)==min(nchar(CDSCHROM))], by = 'GENEID')
txgenes[, .SD[.N>1], by = 'GENEID']
setdiff(keys1, txgenes$GENEID)

# Map chr using ensdb
keytypes(ensdb)
columns(ensdb)
ensgenes <- select(ensdb, keys = keys1, keytype='ENTREZID', columns = c('SEQNAME'))
ensgenes %<>% as.data.table()
ensgenes %<>% extract(, .SD[nchar(SEQNAME)==min(nchar(SEQNAME))], by = 'ENTREZID')
ensgenes[, .SD[length(unique(SEQNAME))>1], by = 'ENTREZID']
ensgenes %<>% extract(!(ENTREZID==266 & SEQNAME=='Y'))
ensgenes[, .SD[length(unique(SEQNAME))>1], by = 'ENTREZID']
ensgenes[, table(SEQNAME)]
ensgenes[, SEQNAME := paste0('chr', SEQNAME)]
ensgenes[, .SD[.N>1], by = 'ENTREZID']
setdiff(keys1, ensgenes$ENTREZID)

keys2 <- unique(as.character(brunello$`Target Gene Symbol`))
ensgenenames <- select(ensdb, keys = keys2, keytype='SYMBOL', columns = c('SEQNAME'))
setdiff(keys2, ensgenenames$SYMBOL)




brunello %>% merge(tmp, by.x = 'Target Gene ID', by.y = 'ENTREZID', sort = FALSE, all.x = TRUE)


tmp %<>% as.data.table()
tmp %<>% extract(!is.na(CDSCHROM))
tmp %<>% extract(, .SD[nchar(CDSCHROM) == min(nchar(CDSCHROM))], by = 'GENEID')
tmp %>% setnames(c('CDSCHROM', 'CDSSTRAND'), c('chrom', 'strand'))
tmp[, .N, by = 'GENEID'][N>1]

brunello %>% merge(tmp, by.x = 'Target Gene ID', by.y = 'GENEID', sort = FALSE, all.x = TRUE)

exons

tx

columns(tx)
columns(ensdb)
select(tx, keys='1', columns = c("TXCHROM", "TXEND", "TXID", "TXNAME", "TXSTART", "TXSTRAND", "TXTYPE"),              keytype = 'GENEID')
select(tx, keys='1', columns = c("TXID", "EXONCHROM","EXONEND", "EXONID", "EXONNAME", "EXONRANK", "EXONSTART", "EXONSTRAND"), keytype = 'GENEID')
select(tx, keys='1', columns = c("TXCHROM", "TXEND", "TXID", "TXNAME", "TXSTART", "TXSTRAND", "TXTYPE"),              keytype = 'GENEID')


ensdb
exons <- ensembldb::exons(ensdb)
txexons <- exons(tx) # 687K
brunello$`Exon Number`
txexons$exon_id
exons # 830K

bs <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
gr <- char_to_granges(c(A1BG= 'chr19:58351502:-'), bs)
gr %<>% extend(-16, +3)
getSeq(bs, gr)
brunello[1]$`sgRNA Target Sequence`

