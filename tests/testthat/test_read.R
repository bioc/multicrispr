

context('get_bsgenome')
test_that('get_bsgenome returns a BSgenome', {
    bedfile  <- system.file('extdata/SRF.bed', package='multicrispr')
    granges  <- read_bed(bedfile, 'mm10')
    expect_s4_class(get_bsgenome(granges), 'BSgenome')
})


context('read_bed')
test_that('read_bed returns a GRanges', {
    bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
    expect_s4_class(read_bed(bedfile, 'mm10', plot = FALSE), 'GRanges')
})
