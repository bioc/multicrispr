context('granges <-> data.table')
test_that('Conversion granges <-> data.table works', {
    
    dt <- data.table::data.table(
        seqnames = 'chr1',
        start    = c( 1,  5, 5),
        end      = c(10, 15, 15),
        strand   = c('+', '+', '-')
    )
    
    bsgenome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
    expect_s4_class(as.granges(dt, bsgenome), 'GRanges')
})

