

context('get_bsgenome')
test_that('get_bsgenome returns a BSgenome', {
    bedfile  <- system.file('extdata/SRF.bed', package='multicrispr')
    bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
    granges  <- read_bed(bedfile, bsgenome)
    expect_s4_class(get_bsgenome(granges), 'BSgenome')
})


context('read_bed')
test_that('read_bed returns a GRanges', {
    bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
    bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
    expect_s4_class(read_bed(bedfile, bsgenome, plot = FALSE), 'GRanges')
})
