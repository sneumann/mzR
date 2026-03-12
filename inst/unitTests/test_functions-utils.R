test_validateHeader <- function() {
    f <- MsDataHub::PestMix1_DDA.mzML()
    o <- openMSfile(f, backend = "pwiz")    
    hdr <- header(o)
    mzR::close(o)
    checkTrue(is(mzR:::.validateHeader(hdr), "data.frame"))
    hdr_2 <- mzR:::.validateHeader(hdr)
    checkTrue(is.character(hdr_2$spectrumId))
    hdr_2 <- mzR:::.validateHeader(hdr[, colnames(hdr) != "spectrumId"])
    checkTrue(is.character(hdr_2$spectrumId))
    checkEquals(hdr_2$spectrumId, paste0("scan=", hdr_2$acquisitionNum))    
    ## Check errors.
    res <- mzR:::.validateHeader(hdr[, 1:5])
    checkTrue(is.character(res))
    hdr2 <- hdr
    hdr2[, 3] <- as.character(hdr2[, 3])
    res <- mzR:::.validateHeader(hdr2)
    checkTrue(is.character(res))
    res <- mzR:::.validateHeader(4)
    checkTrue(is.character(res))
}

test_validSpectrumList <- function() {
    f <- MsDataHub::PestMix1_DDA.mzML()
    o <- openMSfile(f, backend = "pwiz")
    pks <- peaks(o)
    mzR::close(o)
    checkTrue(mzR:::.validSpectrumList(pks))
    ## Check errors.
    res <- mzR:::.validSpectrumList(4)
    checkTrue(is.character(res))
    res <- mzR:::.validSpectrumList(list(4, 2, 4))
    checkTrue(is.character(res))
    pks[[length(pks)]] <- 4
    res <- mzR:::.validSpectrumList(pks)
    checkTrue(is.character(res))
}

test_check_software_processing <- function() {
    checkException(mzR:::.check_software_processing("a"))
    res <- mzR:::.check_software_processing(c("mzR", "1.0.0", "MS:-1"))
    checkEquals(class(res), "list")
    checkEquals(res, list(c("mzR", "1.0.0", "MS:-1")))
    checkException(mzR:::.check_software_processing(c("mzR", "1.0.0")))
    checkException(mzR:::.check_software_processing(c(3, 5)))
}
