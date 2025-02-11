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
    chr <- chromatogram(x, 138L)
    checkTrue(is.data.frame(chr))
    checkIdentical(colnames(chr), c("rtime", "intensity"))
    chr1 <- chromatogram(x, 138L, drop = FALSE)
    checkTrue(is.list(chr1))
    checkEquals(chr, chr1[[1L]])
    close(x)
}

test_chromatograms2 <- function() {
    f <- proteomics(full.names = TRUE,
                    pattern = "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz")
    x <- openMSfile(f, backend = "pwiz")
    checkIdentical(nChrom(x), 1L)
    checkIdentical(tic(x), chromatogram(x, 1L))
    checkIdentical(nrow(tic(x)), 7534L)
    close(x)
}

test_chromatogramHeader_indexing <- function() {
    f <- proteomics(full.names = TRUE, pattern = "MRM")
    x <- openMSfile(f, backend = "pwiz")
    tic <- chromatogramHeader(x, 1)
    tic1 <- chromatogramHeader(x, 1:1)
    # next chromatogram after TIC (if we assume TIC is always at index 1)
    ch2 <- chromatogramHeader(x, 2)
    checkEquals(colnames(ch2), c("chromatogramId", "chromatogramIndex", "polarity",
                                "precursorIsolationWindowTargetMZ",
                                "precursorIsolationWindowLowerOffset",
                                "precursorIsolationWindowUpperOffset",
                                "precursorCollisionEnergy",
                                "productIsolationWindowTargetMZ",
                                "productIsolationWindowLowerOffset",
                                "productIsolationWindowUpperOffset"))
    checkEquals(tic$chromatogramId, "TIC")
    checkEquals(tic1$chromatogramId, "TIC")
    checkTrue(!is.na(ch2$precursorIsolationWindowTargetMZ))
    checkTrue(!is.na(ch2$productIsolationWindowTargetMZ))
    checkEquals(nrow(ch2), 1)
    # test index below lower bound should raise error
    tryCatch({
        chromatogramHeader(x, 0)
        checkTrue(FALSE)
    }, error = function(e) {
        checkTrue(grepl("Index out of bound", e$message))
    })
    # check last chromatogram in bounds
    ch138 <- chromatogramHeader(x, 138)
    checkTrue(startsWith(ch138$chromatogramId, "- SRM SIC Q1=808.1 Q3=407.996"))
    # test index above upper bound should raise error
    tryCatch({
        chromatogramHeader(x, 139)
        checkTrue(FALSE)
    }, error = function(e) {
        checkTrue(grepl("Index out of bound", e$message))
    })

    all_chrom <- chromatogram(x)
    chrom2 <- chromatogram(x, 2)
    checkEquals(all_chrom[[2]], chrom2)

    close(x)
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
