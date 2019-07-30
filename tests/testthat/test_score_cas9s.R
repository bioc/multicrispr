skip_if_no_azimuth <- function(){
    have_azimuth <- reticulate::py_module_available('azimuth')
    if (!have_azimuth) skip("azimuth not available for testing")
}

testthat("score_cas9s works as expected", {
    skip_if_no_azimuth()
})