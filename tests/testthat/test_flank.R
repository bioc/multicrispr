
# Create example granges
bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
granges <-  GenomicRanges::GRanges(
                'chr1', '100-200', strand = '-', targetname = 'T01',
                seqinfo = seqinfo(bsgenome))
chrlength <- GenomeInfoDb::seqlengths(bsgenome)[['chr1']]

# Test
context('up_flank and down_flank')

test_that('up_flank works', {
    expect_equal(   
        start( up_flank(granges, -10, -1)),
        end(granges) + 1 )
    
    expect_equal(
        end( up_flank(granges, -10, -1)),
        end(granges) + 10 )
})

test_that('down_flank works', {
    expect_equal(
        start( down_flank(granges, 1, 10)),
        start(granges) - 10)
    
    expect_equal(
        end( down_flank(granges, 1, 10)),
        start(granges) - 1)
})

test_that('up_flank throws error for coordinates < 1', {
    expect_error(
        up_flank(granges,  +500))
})

test_that('down_flank warns error for coordinates < 1', {
    expect_warning(
        down_flank(granges, -500))
})

test_that('up_flank warns for coordinates > chrlength', {
    expect_warning(
        up_flank( granges, 1, 1 + chrlength))
})

test_that('down_flank warns for coordinates > chrlength', {
    expect_warning(
        down_flank(granges, 1, 1 + chrlength))
})

