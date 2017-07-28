## Testing to write stuff.

test_copyWriteMSData <- function() {
    library(msdata)
    library(mzR)
    library(RUnit)
    test_folder = tempdir()

    ## INPUT: mzXML
    orig_file <- system.file("threonine", "threonine_i2_e35_pH_tree.mzXML",
                        package = "msdata")
    mzxml <- openMSfile(orig_file, backend = "pwiz")
    pks <- peaks(mzxml)
    hdr <- header(mzxml)
    ii <- mzR::instrumentInfo(mzxml)
    mzR::close(mzxml)

    ## OUTPUT: mzML
    fnew <- paste0(test_folder, "test_copyWrite.mzML")
    mzR:::copyWriteMSData(filename = fnew, original_file = orig_file,
                          header = hdr, data = pks, backend = "pwiz")
    ## Check content is same
    mzml_new <- openMSfile(fnew, backend = "pwiz")
    pks_new <- peaks(mzml_new)
    hdr_new <- header(mzml_new)
    ii_new <- mzR::instrumentInfo(mzml_new)
    mzR::close(mzml_new)
    checkEquals(pks_new, pks)
    checkEquals(hdr_new, hdr)
    checkEquals(ii, ii_new)
    
    ## OUTPUT: mzXML
    fnew <- paste0(test_folder, "test_copyWrite.mzXML")
    mzR:::copyWriteMSData(filename = fnew, original_file = orig_file,
                          header = hdr, data = pks, backend = "pwiz",
                          outformat = "mzxml")
    ## Check content is same
    mzml_new <- openMSfile(fnew, backend = "pwiz")
    pks_new <- peaks(mzml_new)
    hdr_new <- header(mzml_new)
    ii_new <- mzR::instrumentInfo(mzml_new)
    mzR::close(mzml_new)
    checkEquals(pks_new, pks)
    checkEquals(hdr_new, hdr)
    checkEquals(ii, ii_new)

    ## Save as mgf
    ## fnew <- paste0(test_folder, "test_copyWrite.mgf")
    ## mzR:::copyWriteMSData(filename = fnew, original_file = orig_file,
    ##                       header = hdr, data = pks, backend = "pwiz",
    ##                       outformat = "mgf")
    ## ## Check content is same
    ## mzml_new <- openMSfile(fnew, backend = "pwiz")
    ## pks_new <- peaks(mzml_new)
    ## hdr_new <- header(mzml_new)
    ## mzR::close(mzml_new)
    ## checkEquals(pks_new, pks)
    ## checkEquals(hdr_new, hdr)
    
    ## Now, let's pick selected spectra instead.
    hdr_sub <- hdr[c(1, 3, 5), ]
    pks_sub <- pks[c(1, 3, 5)]
    fnew <- paste0(test_folder, "test_copyWrite.mzML")
    ## index is not OK after subsetting
    checkException(mzR:::copyWriteMSData(filename = fnew,
                                         original_file = orig_file,
                                         header = hdr_sub, data = pks_sub,
                                         backend = "pwiz"))
    hdr_sub$seqNum <- seq_len(nrow(hdr_sub))
    mzR:::copyWriteMSData(filename = fnew, original_file = orig_file,
                          header = hdr_sub, data = pks_sub, backend = "pwiz")
    ## Check content is same
    mzml_new <- openMSfile(fnew, backend = "pwiz")
    pks_new <- peaks(mzml_new)
    hdr_new <- header(mzml_new)
    mzR::close(mzml_new)
    rownames(hdr_new) <- NULL
    rownames(hdr_sub) <- NULL
    checkEquals(pks_new, pks_sub)
    ## NOTE: percursorScanNum AND acquisitionNum is no longer the same!
    checkException(checkEquals(hdr_new, hdr_sub))

    ## Check errors
    ## wrong header.
    ## wrong spectra.
    ## wrong data processing.
    checkException(mzR:::copyWriteMSData(filename = fnew,
                                         original_file = orig_file,
                                         header = pks, data = hdr,
                                         backend = "pwiz"))
    checkException(mzR:::copyWriteMSData(filename = fnew,
                                         original_file = orig_file,
                                         header = hdr, data = hdr,
                                         backend = "pwiz"))
    checkException(mzR:::copyWriteMSData(filename = fnew,
                                         original_file = orig_file,
                                         header = hdr, data = pks,
                                         backend = "Ramp"))
    checkException(mzR:::copyWriteMSData(filename = fnew,
                                         original_file = "somefile",
                                         header = hdr, data = pks,
                                         backend = "pwiz"))
    checkException(mzR:::copyWriteMSData(filename = fnew,
                                         original_file = orig_file,
                                         header = hdr, data = pks,
                                         backend = "pwiz",
                                         software_processing = c("other")))
    
    ## INPUT: mzML
    orig_file <- system.file("proteomics",
                             "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz",
                             package = "msdata")
    fl <- openMSfile(orig_file, backend = "pwiz")
    pks <- peaks(fl)
    hdr <- header(fl)
    ii <- mzR::instrumentInfo(fl)
    mzR::close(fl)

    ## OUTPUT: mzML
    fnew <- paste0(test_folder, "test_copyWrite.mzML")
    mzR:::copyWriteMSData(filename = fnew, original_file = orig_file,
                          header = hdr, data = pks, backend = "pwiz")
    ## Check content is same
    mzml_new <- openMSfile(fnew, backend = "pwiz")
    pks_new <- peaks(mzml_new)
    hdr_new <- header(mzml_new)
    ii_new <- mzR::instrumentInfo(mzml_new)
    mzR::close(mzml_new)
    checkEquals(pks_new, pks)
    ## acquisitionNum and precursorScanNum will be different, replace with
    ## factors - order and all has to be the same though.
    hdr$acquisitionNum <- as.integer(factor(hdr$acquisitionNum))
    hdr$precursorScanNum <- as.integer(factor(hdr$precursorScanNum))
    hdr_new$acquisitionNum <- as.integer(factor(hdr_new$acquisitionNum))
    hdr_new$precursorScanNum <- as.integer(factor(hdr_new$precursorScanNum))
    checkEquals(hdr_new, hdr)  ## polarity is OK here
    checkEquals(ii, ii_new)
    
    ## OUTPUT: mzXML
    fnew <- paste0(test_folder, "test_copyWrite.mzXML")
    mzR:::copyWriteMSData(filename = fnew, original_file = orig_file,
                          header = hdr, data = pks, backend = "pwiz",
                          outformat = "mzxml")
    ## Check content is same
    mzml_new <- openMSfile(fnew, backend = "pwiz")
    pks_new <- peaks(mzml_new)
    hdr_new <- header(mzml_new)
    ii_new <- mzR::instrumentInfo(mzml_new)
    mzR::close(mzml_new)
    checkEquals(pks_new, pks)
    ## acquisitionNum and precursorScanNum will be different, replace with
    ## factors - order and all has to be the same though.
    hdr$acquisitionNum <- as.integer(factor(hdr$acquisitionNum))
    hdr$precursorScanNum <- as.integer(factor(hdr$precursorScanNum))
    hdr_new$acquisitionNum <- as.integer(factor(hdr_new$acquisitionNum))
    hdr_new$precursorScanNum <- as.integer(factor(hdr_new$precursorScanNum))
    rt_col <- which(colnames(hdr) == "retentionTime")
    checkEquals(hdr[, rt_col], hdr_new[, rt_col], tolerance = 0.01)
    hdr$injectionTime <- 0  ## injectionTime export not supported.
    checkEquals(hdr[ , -rt_col], hdr_new[ , -rt_col])
    ## checkEquals(ii, ii_new)

    ## Other mzML:
    test_file <- system.file("microtofq", "MM14.mzML", package = "msdata")
    in_file <- openMSfile(test_file, backend = "pwiz")
    hdr <- header(in_file)
    pks <- peaks(in_file)
    ii <- mzR::instrumentInfo(in_file)
    mzR::close(in_file)
    
    ## mzML
    out_file <- paste0(test_folder, "test_copyWrite.mzML")
    mzR:::copyWriteMSData(filename = out_file, original_file = test_file,
                          header = hdr, data = pks,
                          software_processing = c("MSnbase", "2.3.8"))
    in_file <- openMSfile(out_file, backend = "pwiz")
    hdr_2 <- header(in_file)
    pks_2 <- peaks(in_file)
    ii_2 <- mzR::instrumentInfo(in_file)
    mzR::close(in_file)
    checkEquals(hdr, hdr_2)
    checkEquals(pks, pks_2)
    checkEquals(ii, ii_2)
    
    ## mzXML output:
    out_file <- paste0(test_folder, "test_copyWrite.mzXML")
    mzR:::copyWriteMSData(filename = out_file, original_file = test_file,
                          header = hdr, data = pks, outformat = "mzXML",
                          software_processing = c("MSnbase", "2.3.8"))
    in_file <- openMSfile(out_file, backend = "pwiz")
    hdr_2 <- header(in_file)
    pks_2 <- peaks(in_file)
    mzR::close(in_file)
    checkEquals(hdr, hdr_2)
    checkEquals(pks, pks_2)
    checkEquals(ii, ii_2)
}

test_writeMSData <- function() {
    library(msdata)
    library(mzR)
    library(RUnit)
    test_folder = tempdir()
    ## Input: mzXML
    test_file <- system.file("threonine", "threonine_i2_e35_pH_tree.mzXML",
                             package = "msdata")
    in_file <- openMSfile(test_file, backend = "pwiz")
    hdr <- header(in_file)
    pks <- peaks(in_file)
    mzR::close(in_file)

    ## mzML
    out_file <- paste0(test_folder, "test_write.mzML")
    mzR:::writeMSData(filename = out_file, header = hdr, data = pks)
    in_file <- openMSfile(out_file, backend = "pwiz")
    hdr_2 <- header(in_file)
    pks_2 <- peaks(in_file)
    mzR::close(in_file)
    checkEquals(hdr, hdr_2)
    checkEquals(pks, pks_2)

    ## Test subsetting.
    hdr_sub <- hdr[c(1, 3, 5), ]
    hdr_sub$seqNum <- 1:nrow(hdr_sub)
    pks_sub <- pks[c(1, 3, 5)]
    mzR:::writeMSData(out_file, header = hdr_sub, data = pks_sub)
    in_file <- openMSfile(out_file)
    hdr_sub_2 <- header(in_file)
    pks_sub_2 <- peaks(in_file)
    mzR::close(in_file)
    checkEquals(hdr, hdr_2)
    checkEquals(pks, pks_2)
    
    ## mzXML output:
    out_file <- paste0(test_folder, "test_write.mzXML")
    mzR:::writeMSData(filename = out_file, header = hdr, data = pks,
                      outformat = "mzXML")
    in_file <- openMSfile(out_file, backend = "pwiz")
    hdr_2 <- header(in_file)
    pks_2 <- peaks(in_file)
    mzR::close(in_file)
    checkEquals(hdr, hdr_2)
    checkEquals(pks, pks_2)

    ## mgf output:
    ## out_file <- paste0(test_folder, "test_write.mgf")
    ## mzR:::writeMSData(filename = out_file, header = hdr, data = pks,
    ##                   outformat = "mgf")
    ## in_file <- openMSfile(out_file, backend = "pwiz")
    ## hdr_2 <- header(in_file)
    ## pks_2 <- peaks(in_file)
    ## mzR::close(in_file)
    ## checkEquals(hdr, hdr_2)
    ## checkEquals(pks, pks_2)

    ## Input: mzML
    test_file <- system.file("proteomics",
                             "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz",
                             package = "msdata")
    in_file <- openMSfile(test_file, backend = "pwiz")
    hdr <- header(in_file)
    pks <- peaks(in_file)
    mzR::close(in_file)
    
    ## mzML
    out_file <- paste0(test_folder, "test_write.mzML")
    mzR:::writeMSData(filename = out_file, header = hdr, data = pks)
    in_file <- openMSfile(out_file, backend = "pwiz")
    hdr_2 <- header(in_file)
    pks_2 <- peaks(in_file)
    mzR::close(in_file)
    checkEquals(hdr, hdr_2)
    checkEquals(pks, pks_2)

    ## mzXML output:
    out_file <- paste0(test_folder, "test_write.mzXML")
    mzR:::writeMSData(filename = out_file, header = hdr, data = pks,
                      outformat = "mzXML")
    in_file <- openMSfile(out_file, backend = "pwiz")
    hdr_2 <- header(in_file)
    pks_2 <- peaks(in_file)
    mzR::close(in_file)
    rt_col <- which(colnames(hdr) == "retentionTime")
    checkEquals(hdr[, rt_col], hdr_2[, rt_col], tolerance = 0.01)
    hdr$injectionTime <- 0  ## injectionTime export not supported.
    checkEquals(hdr[ , -rt_col], hdr_2[ , -rt_col])
    checkEquals(pks, pks_2)

    ## Other mzML:
    test_file <- system.file("microtofq", "MM14.mzML", package = "msdata")
    in_file <- openMSfile(test_file, backend = "pwiz")
    hdr <- header(in_file)
    pks <- peaks(in_file)
    mzR::close(in_file)
    
    ## mzML
    out_file <- paste0(test_folder, "test_write.mzML")
    mzR:::writeMSData(filename = out_file, header = hdr, data = pks)
    in_file <- openMSfile(out_file, backend = "pwiz")
    hdr_2 <- header(in_file)
    pks_2 <- peaks(in_file)
    mzR::close(in_file)
    checkEquals(hdr, hdr_2)
    checkEquals(pks, pks_2)

    ## mzXML output:
    out_file <- paste0(test_folder, "test_write.mzXML")
    mzR:::writeMSData(filename = out_file, header = hdr, data = pks,
                      outformat = "mzXML")
    in_file <- openMSfile(out_file, backend = "pwiz")
    hdr_2 <- header(in_file)
    pks_2 <- peaks(in_file)
    mzR::close(in_file)
    checkEquals(hdr, hdr_2)
    checkEquals(pks, pks_2)

    ## mzData:
    test_file <- system.file("iontrap", "extracted.mzData", package = "msdata")
    in_file <- openMSfile(test_file, backend = "Ramp")
    hdr <- header(in_file)
    pks <- peaks(in_file)
    mzR::close(in_file)

    ## mzML
    out_file <- paste0(test_folder, "test_write.mzML")
    mzR:::writeMSData(filename = out_file, header = hdr, data = pks)
    in_file <- openMSfile(out_file, backend = "pwiz")
    hdr_2 <- header(in_file)
    pks_2 <- peaks(in_file)
    mzR::close(in_file)
    checkEquals(hdr, hdr_2)
    checkEquals(pks, pks_2)

    ## mzXML output:
    out_file <- paste0(test_folder, "test_write.mzXML")
    mzR:::writeMSData(filename = out_file, header = hdr, data = pks,
                      outformat = "mzXML")
    in_file <- openMSfile(out_file, backend = "pwiz")
    hdr_2 <- header(in_file)
    pks_2 <- peaks(in_file)
    mzR::close(in_file)
    checkEquals(hdr, hdr_2)
    checkEquals(pks, pks_2)
}
