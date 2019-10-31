
# Create example granges
bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
granges <-  GenomicRanges::GRanges('chr1', '100-200', strand = '-', 
                    seqinfo = BSgenome::seqinfo(bsgenome))
chrlength <- GenomeInfoDb::seqlengths(bsgenome)[['chr1']]

# Test
context('left_flank and right_flank')

test_that('left_flank works', {
    expect_equal(   
        GenomicRanges::start( left_flank(granges, -10, -1, plot = FALSE)),
        GenomicRanges::start(granges) - 10 )
    
    expect_equal(
        GenomicRanges::end( left_flank(granges, -10, -1, plot = FALSE)),
        GenomicRanges::start(granges) - 1 )
})

test_that('right_flank works', {
    expect_equal(
        GenomicRanges::start( right_flank(granges, 1, 10, plot = FALSE)),
        GenomicRanges::end(granges) + 1)
    
    expect_equal(
        GenomicRanges::end( right_flank(granges, 1, 10, plot = FALSE)),
        GenomicRanges::end(granges) + 10)
})

test_that('left_flank warns for coordinates < 1', {
    expect_warning(
        left_flank(granges,  -500, plot = FALSE))
})

test_that('right_flank warns for coordinates < 1', {
    expect_warning(
        right_flank(granges, -500, plot = FALSE))
})

test_that('left_flank warns for coordinates > chrlength', {
    expect_warning(
        left_flank( granges, 1, 1 + chrlength, plot = FALSE))
})

test_that('right_flank warns for coordinates > chrlength', {
    expect_warning(
        right_flank(granges, 1, 1 + chrlength, plot = FALSE))
})

