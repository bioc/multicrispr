context('get_bsgenome()')
test_that('get_bsgenome returns a BSgenome', {
    bedfile  <- system.file('extdata/SRF.bed', package='multicrispr')
    bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
    granges  <- read_bed(bedfile, bsgenome)
    expect_s4_class(get_bsgenome(granges), 'BSgenome')
})




context('as.granges()')
test_that('as.granges returns a GenomicRanges::GRanges', {
    bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
    dt  <-  data.table::data.table( seqnames = 'chr1',
                                    start    = c( 1,  5, 5),
                                    end      = c(10, 15, 15),
                                    strand   = c('+', '+', '-'))
    expect_s4_class(as.granges(dt, bsgenome), 'GRanges')
})




# context('read_bed')
# test_that('read_bed returns a data.table', {
#     bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
#     bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
#     read_bed(bedfile, bsgenome)
# })
