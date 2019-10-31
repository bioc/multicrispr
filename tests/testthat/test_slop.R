# Create example granges
bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
granges <-  GenomicRanges::GRanges('chr1', 
                    '100-200',
                    strand   = '-',
                    seqinfo  = BSgenome::seqinfo(bsgenome))
chrlength <- GenomeInfoDb::seqlengths(bsgenome)[['chr1']]

# Test
context('extend')

test_that('extend works', {
    expect_equal(
        GenomicRanges::start(extend(granges, -5, 5, plot = FALSE)),  
        GenomicRanges::start(granges) - 5)
    expect_equal(  
        GenomicRanges::end(extend(granges, -5, 5, plot = FALSE)),
        GenomicRanges::end(granges) + 5)
})

test_that('extend warns for coordinates < 1', {
    expect_warning(extend(granges,  -500))
})

test_that('extend warns for coordinates > chrlength', {
    expect_warning(extend(granges, 1, 1 + chrlength, plot = FALSE))
})
