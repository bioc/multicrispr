

context('bed_to_granges')
test_that('bed_to_granges returns a GRanges', {
    bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
    expect_s4_class(bed_to_granges(bedfile, 'mm10', plot = FALSE), 'GRanges')
})
