# Load functions
require(reticulate)
use_condaenv('azienv')
import('azimuth')
require(AnnotationDbi)
require(BSgenome)
require(multicrispr)
require(magrittr)
require(biomaRt)
require(data.table)
require(stringi)

# Load Brunello library
brunello_url <- 'https://www.addgene.org/static/cms/filer_public/8b/4c/8b4c89d9-eac1-44b2-bb2f-8fea95672705/broadgpp-brunello-library-contents.txt'
brunello_file <- '../multicrisprout/brunello/brunello.txt'
#download.file(brunello_url, brunello_file)
brunello <- data.table::fread(brunello_file)
brunello %<>% extract(!is.na(`Target Gene ID`)) # rm pos controls
brunello[, .(N = length(unique(`Target Transcript`))), by = 'Target Gene ID'][, table(N)]    #   1 transcr per gene
brunello[, .(N = length(unique(`Exon Number`))),       by = 'Target Gene ID'][, table(N)]    #   1-4 exons per gene
brunello[, .(N = length(unique(`Exon Number`))),       by = 'Target Transcript'][, table(N)] #             per transcript
brunello[`Target Gene ID`==1]
brunello$`Target Transcript` %<>% substr(1, nchar(.)-2)

# Add strand
ensmart <- biomaRt::useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
strandt <- getBM(attributes = c('refseq_mrna', 'strand'), filters = 'refseq_mrna', values = unique(brunello$`Target Transcript`), mart = ensmart)
strandt %<>% as.data.table()
strandt %<>% extract(, .SD[.N==1], by = 'refseq_mrna') # rm the ones that don't map uniquely to a strand
brunello %<>% merge.data.table(strandt, by.x = 'Target Transcript', by.y = 'refseq_mrna') # 76 441 -> 75 260

# Add seqnames
brunello$`Target Gene ID`
chromdt <- getBM(attributes = c('refseq_mrna', 'chromosome_name'), filters = 'refseq_mrna', values = unique(brunello$`Target Transcript`), mart = ensmart)
chromdt %<>% as.data.table()
chromdt %<>% extract(, .SD[nchar(chromosome_name) == min(nchar(chromosome_name))], by = 'refseq_mrna')
chromdt[, table(chromosome_name)]
chromdt %<>% extract(, .SD[.N>1], by = 'chromosome_name')
chromdt
brunello %<>% merge.data.table(chromdt, by.x = 'Target Transcript', by.y = 'refseq_mrna') # 76 441 -> 75 260 -> 75 244

# Find spacer ranges
brunello %>% setnames(c('Strand', 'strand', 'chromosome_name'), c('sense', 'numericstrand', 'seqnames'))
brunello[, strand := vapply(numericstrand, multicrispr:::csign, character(1)) ]
brunello[numericstrand==+1 & brunello$sense=='antisense', strand := '-']
brunello[numericstrand==-1 & brunello$sense=='antisense', strand := '+']

brunello %>% setnames(c('Position of Base After Cut (1-based)'), c( 'start'))
brunello %>% extract(, end := start)
brunello$cutsite <- brunello$start
bsgenome <- BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38

br <- brunello %>% dt2gr(seqinfo(bsgenome)[standardChromosomes(bsgenome)]) %>% as('GRanges')
br[strand(br)=='+'] %<>% extend(-17, +2)
br[strand(br)=='-'] %<>% extend(-16, +3)
br$crisprspacer <- BSgenome::getSeq(bsgenome, br, as.character=TRUE)
br %<>% subset(`sgRNA Target Sequence`==crisprspacer) # Drop 12 non-matches (due to sequence changes). 75 244 -> 75 232

# Are context seqs identical? Yes!
br$context <- br %>% extend(-4,+6) %>% getSeq(bsgenome, ., as.character=TRUE)
subset(br, context == `Target Context Sequence`)

# find_spacers
br$`Target Gene ID` <- br$`Genomic Sequence` <- br$sense <- br$`sgRNA Target Sequence` <- br$`Target Context Sequence` <- NULL
br$numericstrand <- br$cutsite <- br$`Exon Number` <- br$context <- NULL
names(mcols(br)) %<>% stringi::stri_replace_first_fixed('PAM Sequence', 'crisprpam')
names(mcols(br)) %<>% stringi::stri_replace_first_fixed('Rule Set 2 score', 'Doench2016')

br %<>% sort(ignore.strand=TRUE)
br$targetname <- names(br) <- sprintf('br%05d', seq_along(br))

spacers <- br %>% extend(0, +3) %>% multicrispr::find_spacers(bsgenome, plot = FALSE)

    # All brunello spacers are found
    # 4 813 additional spacers are found on rev strands
    spcoords <- unname(as.character(granges(spacers)))
    brcoords <- unname(as.character(granges(br)))
    setdiff(brcoords, spcoords) %>% length()
    setdiff(spcoords, brcoords) %>% length()
    autonomics.plot::plot_venn(list(multicrispr = spcoords, brunello = brcoords))
    brunello_spacers <- spacers[spcoords %in% brcoords]
    names(mcols(brunello_spacers)) %<>% stringi::stri_replace_first_fixed('Doench2016', 'BrunelloDoench2016')
    
    brunello_spacers[1:10] %>% add_efficiency(bsgenome, 'Doench2016', plot = FALSE)
    brunello_spacers %<>% add_efficiency(bsgenome, 'Doench2016', plot = FALSE)
    plot(brunello_spacers$BrunelloDoench2016, brunello_spacers$Doench2016)
    plot(density(brunello_spacers$Doench2016))


# Find overlapping exons
ensdb <- multicrispr::EnsDb.Hsapiens.v98()
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
exons1 <- ensembldb::exons(ensdb)
#exons1 <- exons(txdb)
chromosomes <- setdiff(standardChromosomes(exons1), 'MT')
exons1 %<>% subset(seqnames(exons1) %in% chromosomes)
exons1 %<>% gr2dt() %>% dt2gr(seqinfo(exons1)[chromosomes])
br %<>% gr2dt() %>% dt2gr(seqinfo(br)[chromosomes])

res <- findOverlaps(exons1, extend(br, 0, +3), ignore.strand = TRUE, minoverlap = 23)
exons1 %<>% extract(queryHits(res))
exons1$targetname <- names(br)[subjectHits(res)]
length(unique(exons1$targetname))
br %<>% extract(unique(subjectHits(res)))
length(br)

exonsdt <- gr2dt(exons1) 
exonsdt[, length(unique(targetname))]
exonsdt %<>% extract(, .SD[width == min(width)], by = 'targetname')  # 75 229 / 
exonsdt %<>% extract(, .SD[1], by = 'targetname')

exons1 <- exonsdt %>% dt2gr(seqinfo(exons1))

all(start(exons1) <= start(extend(br, 0, +3)))
all(  end(exons1) >=   end(extend(br, 0, +3)))
exons1
br

exonspacers <- find_spacers(exons1, bsgenome, plot = FALSE) # 2 552 511
brunellocoords   <- unname(as.character(granges(br)))
exonspacercoords <- unname(as.character(granges(exonspacers)))
autonomics.plot::plot_venn(list(brunello = brunellocoords, multicrispr = exonspacercoords))

exonspacers$in_brunello <- exonspacercoords %in% brunellocoords
exonspacers %>% subset(in_brunello==TRUE)
exonspacers %<>% add_context(bsgenome)

exonspacers %<>% add_efficiency(bsgenome, 'Doench2016')

contextseqs <- exonspacers$crisprcontext[1:1000]

score_doench2016_parallel <- function(contextseqs, chunksize = 100){
    
    nchunks <- ceiling(length(contextseqs) / chunksize)
    scores<-BiocParallel::bplapply(
                seq_len(nchunks), 
                function(i){
                   chunkstart <- (i-1)*chunksize + 1
                   chunkend   <- min(i*chunksize, length(contextseqs))
                   multicrispr:::doench2016(contextseqs[chunkstart:chunkend])
                }) %>% 
            unlist()
}



parallel::mclapply()
exonspacers$Doench2016 <- 0
exonspacers$Doench2016[     1:100000] <- multicrispr:::doench2016(exonspacers$crisprcontext[     1:100000])
exonspacers$Doench2016[100001:200000] <- multicrispr:::doench2016(exonspacers$crisprcontext[100001:200000])
tmp <- exonspacers[1:100000]
tmp %<>% add_efficiency(bsgenome, 'Doench2016')


# add_genome_counts
bsgenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
seqlevelsStyle(br) <- 'UCSC'
br %<>% add_genome_counts(bsgenome, outdir = '~/multicrisprout/brunello')
length(unique(br$G0))
table(br$G0)
barplot(table(br$G0), log = "y", xlab = '# Exact genome matches')



brunello #    76 065
spacers  # 4 374 614
spacers %<>% add_genome_counts(  # works (also) on vm!
                mismatches=1, 
                indexedgenomesdir = '../../multicrisprout/indexedgenomes', 
                outdir = '../../multicrisprout')
spacerdt <- as.data.table(spacers)

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

