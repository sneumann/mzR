test_header_all <- function() {
    ## Check that the header call returns the same columns irrespectively of the
    ## backend. Issue #238
    file <- system.file('cdf/ko15.CDF', package = "msdata")
    cdf <- openMSfile(file, backend="netCDF")
    header_cdf <- header(cdf)
    close(cdf)

    file <- system.file("threonine", "threonine_i2_e35_pH_tree.mzXML", package = "msdata")
    mzxml <- openMSfile(file, backend="pwiz")
    header_pwiz <- header(mzxml)
    close(mzxml)

    mzxml <- openMSfile(file, backend="Ramp")
    header_ramp <- header(mzxml)
    close(mzxml)

    checkEquals(colnames(header_cdf), colnames(header_pwiz))
    checkEquals(colnames(header_ramp), colnames(header_pwiz))
}
