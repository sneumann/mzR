test_chromatograms1 <- function() {
    f <- proteomics(full.names = TRUE, pattern = "MRM")
    x <- openMSfile(f, backend = "pwiz")
    checkIdentical(nChrom(x), 138L)
    checkIdentical(tic(x), chromatogram(x, 1L))
    checkIdentical(chromatogram(x), chromatograms(x))
    checkIdentical(nrow(chromatogram(x, 1L)), 85799L)
    checkIdentical(nrow(chromatogram(x, 110L)), 591L)
    checkIdentical(nrow(chromatogram(x, 111L)), 1004L)
    checkIdentical(nrow(chromatogram(x, 112L)), 1004L)
    checkIdentical(nrow(chromatogram(x, 136L)), 527L)
    checkIdentical(nrow(chromatogram(x, 137L)), 567L)
    checkIdentical(nrow(chromatogram(x, 138L)), 567L)
}

test_chromatograms2 <- function() {
    f <- proteomics(full.names = TRUE,
                    pattern = "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz")
    x <- openMSfile(f, backend = "pwiz")
    checkIdentical(nChrom(x), 1L)   
    checkIdentical(tic(x), chromatogram(x, 1L))
    checkIdentical(nrow(tic(x)), 7534L)
}

test_chromatogramHeader <- function() {
    library(mzR)
    library(RUnit)
    library(msdata)
    
    f <- proteomics(full.names = TRUE, pattern = "MRM")
    x <- openMSfile(f)

    chrs <- chromatogram(x)
    ch <- chromatogramHeader(x)
    checkEquals(colnames(ch), c("chromatogramId", "chromatogramIndex", "polarity",
                                "precursorIsolationWindowTargetMZ",
                                "precursorIsolationWindowLowerOffset",
                                "precursorIsolationWindowUpperOffset",
                                "precursorCollisionEnergy",
                                "productIsolationWindowTargetMZ",
                                "productIsolationWindowLowerOffset",
                                "productIsolationWindowUpperOffset"))
    checkEquals(ch$chromatogramId[1], "TIC")
    checkEquals(sum(is.na(ch$precursorIsolationWindowTargetMZ)), 1)
    checkEquals(sum(is.na(ch$productIsolationWindowTargetMZ)), 1)    
    checkEquals(length(chrs), nrow(ch))
    close(x)

    ## Should return only the TIC.
    f <- proteomics(full.names = TRUE, pattern = "MS3")
    x <- openMSfile(f[1])
    ch <- chromatogramHeader(x)
    checkEquals(nrow(ch), 1)
    checkEquals(ch$chromatogramId, "TIC")
    close(x)
}
