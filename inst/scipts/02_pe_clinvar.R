#=============================================================================
# Clinvar is available in VCF (highly formalized) and HGVS (less formalized).
#    The VCF format is highly formalized and well suited for programming, 
#    The HGVS format is ambiguous and not well suited for programming
#=============================================================================


# Clinvar in VCF (preferred)
#---------------------------
devtools::load_all()
localfile <- paste0('~/.multicrispr/clinvar.vcf')
if (!file.exists(localfile)){
    utils::download.file(
        url = 'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz', 
        destfile = paste0(localfile, '.gz'))
}
bsgenome <- BSgenome.Hsapiens.NCBI.GRCh38::BSgenome.Hsapiens.NCBI.GRCh38
#vcf <- VariantAnnotation::readVcf('~/.multicrispr/clinvar.vcf')
vr  <- VariantAnnotation::readVcfAsVRanges('~/.multicrispr/clinvar.vcf')
vr %<>% extract(as.character(seqnames(vr)) %in% seqlevels(bsgenome))
seqlevels(vr) <- seqlevelsInUse(vr)
seqinfo(vr) <- seqinfo(bsgenome)[seqlevels(vr)]
vr
show_sample <- function(clnvc, i=1){ 
    smallvr <- vr[vr$CLNVC == clnvc][i]
    mcols(smallvr) <- NULL
    print(smallvr)
    BSgenome::getSeq(bsgenome, smallvr)
}
sort(table(vr$CLNVC))
show_sample('single_nucleotide_variant')  # 443 238
show_sample('Deletion')                   #  31 472
show_sample('Duplication')                #  14 156
show_sample('Microsatellite')             #   8 506
show_sample('Indel')                      #   3 923
show_sample('Insertion')                  #   2 944
show_sample('Inversion')                  #     226
show_sample('copy_number_loss')           #      13
show_sample('Variation')                  #       6
show_sample('copy_number_gain')           #       2

# ClinVar in HGVS (not well suited for programming)
#--------------------------------------------------
localfile <- '~/.multicrispr/variant_summary.txt.gz'
if (!file.exists(localfile)){
    utils::download.file(
        url = paste0(   
                'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/', 
                'variant_summary.txt.gz'), 
        destfile = localfile)
}
clinvar <- data.table::fread(localfile)
clinvar %<>% extract(Assembly=='GRCh38')                    # 523 K
clinvar %<>% extract(OriginSimple == 'germline')            # 481 K
clinvar[, sort(table(ClinicalSignificance))]
clinvar %<>% extract(ClinicalSignificance == 'Pathogenic')  #  64 K
clinvar[,  c('Assembly', 
             'ClinicalSignificance', 'ClinSigSimple', 
             'Origin', 'OriginSimple', 
             'PhenotypeIDS', 'VariationID', 'SubmitterCategories', 'OtherIDs', 
             'Guidelines', 'TestedInGTR', 'LastEvaluated', 
             'RCVaccession', 'HGNC_ID', 'RS# (dbSNP)', 
             'nsv/esv (dbVar)', 'ChromosomeAccession') := NULL ]
clinvar[, sort(table(Type))]
clinvar[, sort(table(ReviewStatus))]
clinvar[, length(unique(GeneID))]     # 3 790 genes
clinvar[, .(nalleles = .N), by = 'GeneID'][, .(freq = .N), by = nalleles][order(freq)]
clinvar %<>% extract(Stop - Start)

firstrecord <- function(dt, itype, igene=1){
    dt %>% 
    extract(Type==names(rev(sort(table(Type))))[itype]) %>% 
    extract(GeneSymbol==unique(GeneSymbol)[igene]) %>% 
    extract(1)
}
getSeq <- function(dt)  dt %>% 
                        extract(, BSgenome::getSeq(
                                    bsgenome, 
                                    paste0('chr', Chromosome), Start, Stop))
clinvar[, sort(table(Type))]
firstrecord(clinvar, 1)                  # 8 335 snvs
firstrecord(clinvar, 1)    %>% getSeq()
firstrecord(clinvar, 2)                  # 1 950 indels
firstrecord(clinvar, 2, 2)
firstrecord(clinvar, 2, 2) %>% getSeq()
firstrecord(clinvar, 3)                  # 1 604 deletions
firstrecord(clinvar, 3)    %>% getSeq()
firstrecord(clinvar, 4, 3)               # 1 503 duplications
firstrecord(clinvar, 4, 3) %>% getSeq()
firstrecord(clinvar, 5)                  #   377 insertions
firstrecord(clinvar, 5)    %>% getSeq()



clinvar[Type==names(rev(sort(table(Type))))[1]][GeneSymbol==unique(GeneSymbol)[1]][1]

clinvar[Type==names(rev(sort(table(Type))))[2]][GeneSymbol==unique(GeneSymbol)[1]][1]

clinvar
clinvar[Type=='deletion']
cli
