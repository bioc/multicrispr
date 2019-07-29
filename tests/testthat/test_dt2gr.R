context('dt2gr')
test_that('dt2gr and gr2dt work', {
    
    dt <- data.table::data.table(
        chr    = 'chr1',
        start  = c( 1,  5, 5),
        end    = c(10, 15, 15),
        strand = c('+', '+', '-')
    )
    
    expect_s4_class(dt2gr(dt), 'GRanges')
    expect_equal(gr2dt(dt2gr(dt)), dt)
    
})

#' dt %>% dt2gr()
#' dt %>% dt2gr() %>% gr2dt()
