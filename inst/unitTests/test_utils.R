test_hasChromatograms <- function() {
    fl <- system.file("proteomics/MRM-standmix-5.mzML.gz", package = "msdata")
    x <- mzR::openMSfile(fl, backend = "pwiz")
    checkTrue(mzR:::.hasChromatograms(x))
    checkTrue(mzR:::.hasChromatograms(fl))
    close(x)
    
    fl <- system.file("cdf/ko15.CDF", package = "msdata")
    x <- openMSfile(fl, backend = "netCDF")        
    suppressWarnings(checkTrue(!mzR:::.hasChromatograms(x)))
    suppressWarnings(checkTrue(!mzR:::.hasChromatograms(fl)))
    close(x)

    fl <- system.file("sciex/20171016_POOL_POS_1_105-134.mzML",
                      package = "msdata")
    x <- mzR::openMSfile(fl, backend = "pwiz")
    checkTrue(!mzR:::.hasChromatograms(x))
    checkTrue(!mzR:::.hasChromatograms(fl))
    close(x)
}
