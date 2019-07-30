azimuth <- NULL
numpy   <- NULL

.onload <- function(libname, pkgname){
    azimuth <<- reticulate::import("azimuth", delay_load = TRUE)
    numpy   <<- reticulate::import("numpy",   delay_load = TRUE)
}