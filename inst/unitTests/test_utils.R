test_hasChromatograms <- function() {
    fl <- MsDataHub::MRM.standmix.5.mzML()
    x <- mzR::openMSfile(fl)
    checkTrue(mzR:::.hasChromatograms(x))
    checkTrue(mzR:::.hasChromatograms(fl))
    close(x)

    fl <- MsDataHub::ko15.CDF()
    x <- openMSfile(fl)
    suppressWarnings(checkTrue(!mzR:::.hasChromatograms(x)))
    suppressWarnings(checkTrue(!mzR:::.hasChromatograms(fl)))
    close(x)

    fl <- MsDataHub::PestMix1_DDA.mzML()
    x <- openMSfile(fl)
    checkTrue(mzR:::.hasChromatograms(x))
    checkTrue(mzR:::.hasChromatograms(fl))
    close(x)
}
