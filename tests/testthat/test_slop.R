# Create example granges
bsgenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
granges <-  GRanges('chr1', '100-200', strand   = '-',
                                    seqinfo  = BSgenome::seqinfo(bsgenome))
chrlength <- seqlengths(bsgenome)[['chr1']]

# Test
context('slop')

test_that('slop works', {
    expect_equal(start(slop(granges, -5, 5)),  start(granges) - 5)
    expect_equal(  end(slop(granges, -5, 5)),    end(granges) + 5)
})

test_that('slop warns for coordinates < 1', {
    expect_warning(slop(granges,  -500))
})

test_that('slop warns for coordinates > chrlength', {
    expect_warning(slop(granges, 1, 1 + chrlength))
})
