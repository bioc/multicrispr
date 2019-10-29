

context('bed_to_granges')
test_that('bed_to_granges returns a GRanges', {
    bedfile <- system.file('extdata/SRF.bed', package = 'multicrispr')
    txdb <- utils::getFromNamespace('TxDb.Mmusculus.UCSC.mm10.knownGene', 
                                    'TxDb.Mmusculus.UCSC.mm10.knownGene')
    expect_s4_class(bed_to_granges(bedfile, txdb, plot = FALSE), 'GRanges')
})
