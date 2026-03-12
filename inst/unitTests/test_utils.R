test_hasChromatograms <- function() {
    x <- mzR::openMSfile(MsDataHub::MRM.standmix.5.mzML())
    checkTrue(mzR:::.hasChromatograms(x))
    checkTrue(mzR:::.hasChromatograms(fl))
    close(x)
    
    x <- openMSfile(MsDataHub::ko15.CDF())
    suppressWarnings(checkTrue(!mzR:::.hasChromatograms(x)))
    suppressWarnings(checkTrue(!mzR:::.hasChromatograms(fl)))
    close(x)

    x <- mzR::openMSfile(MsDataHub::PestMix1_DDA.mzML())
    checkTrue(mzR:::.hasChromatograms(x))
    checkTrue(mzR:::.hasChromatograms(fl))
    close(x)
}
