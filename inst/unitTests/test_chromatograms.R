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
    f <- proteomics(full.names = TRUE, pattern = "^TMT")
    x <- openMSfile(f, backend = "pwiz")
    checkIdentical(nChrom(x), 1L)   
    checkIdentical(tic(x), chromatogram(x, 1L))
    checkIdentical(nrow(tic(x)), 7534L)
}
