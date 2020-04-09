###
###

.pkgname <- "BSgenome.SarsCov2.UCSC.wuhCor1"

.seqnames <- NULL

.circ_seqs <- NULL

.mseqnames <- NULL

.onLoad <- function(libname, pkgname)
{
    if (pkgname != .pkgname)
        stop("package name (", pkgname, ") is not ",
             "the expected name (", .pkgname, ")")
    extdata_dirpath <- system.file("extdata", package=pkgname,
                                   lib.loc=libname, mustWork=TRUE)

    ## Make and export BSgenome object.
    bsgenome <- BSgenome(
        organism="SARS-Cov-2",
        common_name="Wuhan Coronavirus",
        provider="UCSC",
        provider_version="wuhCor1",
        release_date="Jan. 2020",
        release_name="wuhCor1",
        source_url="http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/bigZips",
        seqnames=.seqnames,
        circ_seqs=.circ_seqs,
        mseqnames=.mseqnames,
        seqs_pkgname=pkgname,
        seqs_dirpath=extdata_dirpath
    )

    ns <- asNamespace(pkgname)

    objname <- pkgname
    assign(objname, bsgenome, envir=ns)
    namespaceExport(ns, objname)

    old_objname <- "SarsCov2"
    assign(old_objname, bsgenome, envir=ns)
    namespaceExport(ns, old_objname)
}

