## Testing to write stuff.

test_copyWriteMSData <- function() {
    test_folder = tempdir()

    ## INPUT: mzXML
    orig_file <- system.file("threonine", "threonine_i2_e35_pH_tree.mzXML",
                             package = "msdata")

    mzML_xsd <- XML::xmlTreeParse(system.file("extdata", "mzML1.1.0.xsd",
                                              package = "mzR"),
                                  isSchema = TRUE, useInternal = TRUE)
    mzML_xsd_idx <- XML::xmlTreeParse(system.file("extdata", "mzML1.1.2_idx.xsd",
                                                  package = "mzR"),
                                      isSchema = TRUE, useInternal = TRUE)

    mzxml <- openMSfile(orig_file, backend = "pwiz")
    pks <- peaks(mzxml)
    hdr <- header(mzxml)
    ii <- mzR::instrumentInfo(mzxml)
    mzR::close(mzxml)

    ## OUTPUT: mzML
    fnew <- paste0(test_folder, "test_copyWrite.mzML")
    mzR::copyWriteMSData(file = fnew, original_file = orig_file,
                         header = hdr, object = pks, backend = "pwiz")
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
    mzR::copyWriteMSData(file = fnew, original_file = orig_file,
                         header = hdr, object = pks, backend = "pwiz",
                         outformat = "mzxml")
    ## Check content is same
    mzml_new <- openMSfile(fnew, backend = "pwiz")
    pks_new <- peaks(mzml_new)
    hdr_new <- header(mzml_new)
    ii_new <- mzR::instrumentInfo(mzml_new)
    mzR::close(mzml_new)
    checkEquals(pks_new, pks)
    ## Don't compare IDs since they are different.
    checkEquals(hdr_new[, colnames(hdr_new) != "spectrumId"],
                hdr[, colnames(hdr) != "spectrumId"])
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
    checkException(mzR:::copyWriteMSData(file = fnew,
                                         original_file = orig_file,
                                         header = hdr_sub, object = pks_sub,
                                         backend = "pwiz"))
    hdr_sub$seqNum <- seq_len(nrow(hdr_sub))
    ## mzML
    mzR::copyWriteMSData(file = fnew, original_file = orig_file,
                         header = hdr_sub, object = pks_sub, backend = "pwiz",
                         outformat = "mzml")
    ## Check content is same
    mzml_new <- openMSfile(fnew, backend = "pwiz")
    pks_new <- peaks(mzml_new)
    hdr_new <- header(mzml_new)
    rownames(hdr_new) <- NULL
    rownames(hdr_sub) <- NULL
    checkEquals(pks_new, pks_sub)
    checkEquals(hdr_new, hdr_sub)
    checkEquals(peaks(mzml_new, 2), pks[[3]])
    mzR::close(mzml_new)
    ## mzXML
    mzR::copyWriteMSData(file = fnew, original_file = orig_file,
                         header = hdr_sub, object = pks_sub, backend = "pwiz",
                         outformat = "mzxml")
    ## Check content is same
    mzml_new <- openMSfile(fnew, backend = "pwiz")
    pks_new <- peaks(mzml_new)
    hdr_new <- header(mzml_new)
    mzR::close(mzml_new)
    rownames(hdr_new) <- NULL
    rownames(hdr_sub) <- NULL
    ## acquisitionNum and precursorScanNum are expected to be different, same
    ## as spectrumId
    hdr_new$acquisitionNum <- as.integer(factor(hdr_new$acquisitionNum))
    hdr_sub$acquisitionNum <- as.integer(factor(hdr_sub$acquisitionNum))
    hdr_new$precursorScanNum <- as.integer(factor(hdr_new$precursorScanNum))
    hdr_sub$precursorScanNum <- as.integer(factor(hdr_sub$precursorScanNum))
    hdr_new$spectrumId <- as.integer(factor(hdr_new$spectrumId))
    hdr_sub$spectrumId <- as.integer(factor(hdr_sub$spectrumId))
    checkEquals(pks_new, pks_sub)
    checkEquals(hdr_new, hdr_sub)

    ## Check errors
    ## wrong header.
    ## wrong spectra.
    ## wrong data processing.
    checkException(mzR::copyWriteMSData(file = fnew,
                                        original_file = orig_file,
                                        header = pks, object = hdr,
                                        backend = "pwiz"))
    checkException(mzR::copyWriteMSData(file = fnew,
                                        original_file = orig_file,
                                        header = hdr, object = hdr,
                                        backend = "pwiz"))
    checkException(mzR::copyWriteMSData(file = fnew,
                                        original_file = orig_file,
                                        header = hdr, object = pks,
                                        backend = "Ramp"))
    checkException(mzR::copyWriteMSData(file = fnew,
                                        original_file = "somefile",
                                        header = hdr, object = pks,
                                        backend = "pwiz"))
    checkException(mzR::copyWriteMSData(file = fnew,
                                        original_file = orig_file,
                                        header = hdr, object = pks,
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
    mzR::copyWriteMSData(file = fnew, original_file = orig_file,
                         header = hdr, object = pks, backend = "pwiz")
    ## Check content is same
    mzml_new <- openMSfile(fnew, backend = "pwiz")
    pks_new <- peaks(mzml_new)
    hdr_new <- header(mzml_new)
    ii_new <- mzR::instrumentInfo(mzml_new)
    checkEquals(pks_new, pks)
    checkEquals(hdr_new, hdr)  ## polarity is OK here
    checkEquals(ii, ii_new)

    checkEquals(peaks(mzml_new, 12), pks[[12]])
    mzR::close(mzml_new)
    ## OUTPUT: mzXML
    fnew <- paste0(test_folder, "test_copyWrite.mzXML")
    suppressWarnings(
        mzR::copyWriteMSData(file = fnew, original_file = orig_file,
                             header = hdr, object = pks, backend = "pwiz",
                             outformat = "mzxml")
    )
    ## Check content is same
    mzml_new <- openMSfile(fnew, backend = "pwiz")
    pks_new <- peaks(mzml_new)
    hdr_new <- header(mzml_new)
    ii_new <- mzR::instrumentInfo(mzml_new)
    mzR::close(mzml_new)
    checkEquals(pks_new, pks)
    ## acquisitionNum and precursorScanNum will be different, replace with
    ## factors - order and all has to be the same though.
    hdr_mod <- hdr
    hdr_mod$acquisitionNum <- as.integer(factor(hdr$acquisitionNum))
    hdr_mod$precursorScanNum <- as.integer(factor(hdr$precursorScanNum))
    hdr_new$acquisitionNum <- as.integer(factor(hdr_new$acquisitionNum))
    hdr_new$precursorScanNum <- as.integer(factor(hdr_new$precursorScanNum))
    rt_col <- colnames(hdr_mod) == "retentionTime"
    checkEquals(hdr_mod[, rt_col], hdr_new[, rt_col], tolerance = 0.01)
    cn <- colnames(hdr)[!(colnames(hdr) %in% c("injectionTime", "retentionTime",
                                               "filterString", "spectrumId"))]
    checkEquals(hdr_mod[, cn], hdr_new[, cn])
    ## checkEquals(ii, ii_new)

    ## Subset. These checks ensure that the scan - precursor scan are mapped
    ## correctly!
    idx <- c(`1003` = 1, `1006` = 4, `1019` = 17, `1021` = 19, `1024` = 22,
             `1026` = 24)
    ## 1: no precursor
    ## 2: 1 as precursor, 1003
    ## 3: no precursor
    ## 4: no precursor
    ## 5: 4 as precursor, 1021
    ## 6: 4 as precursor, 1021
    hdr_sub <- hdr[idx, ]
    pks_sub <- pks[idx]
    hdr_sub$seqNum <- 1:nrow(hdr_sub)
    fnew <- paste0(test_folder, "test_copyWrite.mzML")
    mzR::copyWriteMSData(file = fnew, original_file = orig_file,
                         header = hdr_sub, object = pks_sub, backend = "pwiz",
                         outformat = "mzml")
    ## Check content is same
    mzml_new <- openMSfile(fnew, backend = "pwiz")
    pks_new <- peaks(mzml_new)
    hdr_new <- header(mzml_new)
    ii_new <- mzR::instrumentInfo(mzml_new)
    mzR::close(mzml_new)
    checkEquals(pks_new, pks_sub)
    rownames(hdr_sub) <- NULL
    rownames(hdr_new) <- NULL
    checkEquals(hdr_new, hdr_sub)  ## polarity is OK here
    checkEquals(ii, ii_new)

    ## Subset with mzXML
    fnew <- paste0(test_folder, "test_copyWrite.mzXML")
    suppressWarnings(
        mzR::copyWriteMSData(file = fnew, original_file = orig_file,
                             header = hdr_sub, object = pks_sub, backend = "pwiz",
                             outformat = "mzxml")
    )
    ## Check content is same
    mzml_new <- openMSfile(fnew, backend = "pwiz")
    pks_new <- peaks(mzml_new)
    hdr_new <- header(mzml_new)
    ii_new <- mzR::instrumentInfo(mzml_new)
    mzR::close(mzml_new)
    checkEquals(pks_new, pks_sub)
    rownames(hdr_sub) <- NULL
    rownames(hdr_new) <- NULL
    hdr_sub$acquisitionNum <- as.integer(factor(hdr_sub$acquisitionNum))
    hdr_sub$precursorScanNum <- as.integer(factor(hdr_sub$precursorScanNum))
    hdr_new$acquisitionNum <- as.integer(factor(hdr_new$acquisitionNum))
    hdr_new$precursorScanNum <- as.integer(factor(hdr_new$precursorScanNum))
    rt_col <- colnames(hdr_sub) == "retentionTime"
    checkEquals(hdr_sub[, rt_col], hdr_new[, rt_col], tolerance = 0.01)
    cn <- colnames(hdr_sub)[!(colnames(hdr_sub) %in% c("injectionTime",
                                                       "retentionTime",
                                                       "filterString",
                                                       "spectrumId"))]
    checkEquals(hdr_sub[, cn], hdr_new[, cn])    
    
    ## Other mzML:
    test_file <- system.file("microtofq", "MM14.mzML", package = "msdata")
    in_file <- openMSfile(test_file, backend = "pwiz")
    hdr <- header(in_file)
    pks <- peaks(in_file)
    ii <- mzR::instrumentInfo(in_file)
    mzR::close(in_file)
    
    ## mzML
    out_file <- paste0(test_folder, "test_copyWrite.mzML")
    mzR::copyWriteMSData(file = out_file, original_file = test_file,
                         header = hdr, object = pks,
                         software_processing = c("MSnbase", "2.3.8", "MS:-1"))
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
    mzR::copyWriteMSData(file = out_file, original_file = test_file,
                         header = hdr, object = pks, outformat = "mzXML",
                         software_processing = c("MSnbase", "2.3.8", "MS:-1"))
    in_file <- openMSfile(out_file, backend = "pwiz")
    hdr_2 <- header(in_file)
    pks_2 <- peaks(in_file)
    mzR::close(in_file)
    checkEquals(hdr[, colnames(hdr) != "spectrumId"],
                hdr_2[, colnames(hdr_2) != "spectrumId"])
    checkEquals(pks, pks_2)
    checkEquals(ii, ii_2)
}

test_writeMSData <- function() {    
    mzML_xsd <- XML::xmlTreeParse(system.file("extdata", "mzML1.1.0.xsd",
                                              package = "mzR"),
                                  isSchema = TRUE, useInternal = TRUE)
    mzML_xsd_idx <- XML::xmlTreeParse(system.file("extdata", "mzML1.1.2_idx.xsd",
                                                  package = "mzR"),
                                      isSchema = TRUE, useInternal = TRUE)
    mzXML_xsd_idx <- XML::xmlTreeParse(system.file("extdata",
                                                   "mzXML_idx_3.2.xsd.xml",
                                                   package = "mzR"),
                                       isSchema = TRUE, useInternal = TRUE)

    test_folder = tempdir()
    ## Input: mzXML
    test_file <- system.file("threonine", "threonine_i2_e35_pH_tree.mzXML",
                             package = "msdata")
    in_file <- openMSfile(test_file, backend = "pwiz")
    hdr <- header(in_file)
    pks <- peaks(in_file)
    mzR::close(in_file)

    ## mzML
    out_file <- paste0(test_folder, "/test_write.mzML")
    writeMSData(file = out_file, header = hdr, object = pks)
    in_file <- openMSfile(out_file, backend = "pwiz")
    hdr_2 <- header(in_file)
    pks_2 <- peaks(in_file)
    checkEquals(hdr, hdr_2)
    checkEquals(pks, pks_2)
    checkEquals(peaks(in_file, 13), pks[[13]])
    mzR::close(in_file)
    ## validate mzML:
    doc <- XML::xmlInternalTreeParse(out_file)
    res <- XML::xmlSchemaValidate(mzML_xsd_idx, doc)
#    checkEquals(res$status, 0)
    
    ## Test subsetting.
    hdr_sub <- hdr[c(1, 3, 5), ]
    hdr_sub$seqNum <- 1:nrow(hdr_sub)
    pks_sub <- pks[c(1, 3, 5)]
    writeMSData(pks_sub, out_file, header = hdr_sub)
    in_file <- openMSfile(out_file)
    hdr_sub_2 <- header(in_file)
    pks_sub_2 <- peaks(in_file)
    checkEquals(pks_sub, pks_sub_2)
    checkEquals(peaks(in_file, 3), pks[[5]])
    mzR::close(in_file)
    ## mzXML does not support spectrumId, thus acquisitionNum, precursorScanNum
    ## and spectrumId will be different, but their order and mapping has to be
    ## the same.
    hdr_sub$acquisitionNum <- as.integer(factor(hdr_sub$acquisitionNum))
    hdr_sub_2$acquisitionNum <- as.integer(factor(hdr_sub_2$acquisitionNum))
    hdr_sub$precursorScanNum <- as.integer(factor(hdr_sub$precursorScanNum))
    hdr_sub_2$precursorScanNum <- as.integer(factor(hdr_sub_2$precursorScanNum))
    hdr_sub$spectrumId <- as.integer(factor(hdr_sub$spectrumId))
    hdr_sub_2$spectrumId <- as.integer(factor(hdr_sub_2$spectrumId))
    rownames(hdr_sub) <- NULL
    checkEquals(hdr_sub, hdr_sub_2)
    ## validate mzML:
    doc <- XML::xmlInternalTreeParse(out_file)
    res <- XML::xmlSchemaValidate(mzML_xsd_idx, doc)
    checkEquals(res$status, 0)
    
    ## mzXML output:
    out_file <- paste0(test_folder, "/test_write.mzXML")
    writeMSData(file = out_file, header = hdr, object = pks,
                outformat = "mzXML")
    in_file <- openMSfile(out_file, backend = "pwiz")
    hdr_2 <- header(in_file)
    pks_2 <- peaks(in_file)
    mzR::close(in_file)
    checkEquals(pks, pks_2)
    checkEquals(hdr[, colnames(hdr) != "spectrumId"],
                hdr_2[, colnames(hdr_2) != "spectrumId"])
    hdr_sub <- hdr[c(1, 3, 5), ]
    hdr_sub$seqNum <- 1:nrow(hdr_sub)
    pks_sub <- pks[c(1, 3, 5)]
    writeMSData(file = out_file, header = hdr_sub, object = pks_sub,
                outformat = "mzXML")
    in_file <- openMSfile(out_file)
    hdr_sub_2 <- header(in_file)
    pks_sub_2 <- peaks(in_file)
    mzR::close(in_file)
    checkEquals(pks_sub, pks_sub_2)
    ## mzXML does not support spectrumId, thus acquisitionNum, precursorScanNum
    ## and spectrumId will be different, but their order and mapping has to be
    ## the same.
    hdr_sub$acquisitionNum <- as.integer(factor(hdr_sub$acquisitionNum))
    hdr_sub_2$acquisitionNum <- as.integer(factor(hdr_sub_2$acquisitionNum))
    hdr_sub$precursorScanNum <- as.integer(factor(hdr_sub$precursorScanNum))
    hdr_sub_2$precursorScanNum <- as.integer(factor(hdr_sub_2$precursorScanNum))
    hdr_sub$spectrumId <- as.integer(factor(hdr_sub$spectrumId))
    hdr_sub_2$spectrumId <- as.integer(factor(hdr_sub_2$spectrumId))
    rownames(hdr_sub) <- NULL
    checkEquals(hdr_sub, hdr_sub_2)
    
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
    out_file <- paste0(test_folder, "/test_write.mzML")
    writeMSData(pks, file = out_file, header = hdr)
    in_file <- openMSfile(out_file, backend = "pwiz")
    hdr_2 <- header(in_file)
    pks_2 <- peaks(in_file)
    mzR::close(in_file)
    checkEquals(pks, pks_2)
    ## Don't understand exactly why, but here we do have different
    ## acquisitionNum and precursorScanNum while having the same spectrumID.
    hdr_mod <- hdr
    hdr_mod$acquisitionNum <- as.integer(factor(hdr_mod$acquisitionNum))
    hdr_mod$precursorScanNum <- as.integer(factor(hdr_mod$precursorScanNum))
    hdr_2$acquisitionNum <- as.integer(factor(hdr_2$acquisitionNum))
    hdr_2$precursorScanNum <- as.integer(factor(hdr_2$precursorScanNum))
    checkEquals(hdr_mod, hdr_2)
    ## validate mzML:
    doc <- XML::xmlInternalTreeParse(out_file)
    res <- XML::xmlSchemaValidate(mzML_xsd_idx, doc)
    checkEquals(res$status, 0)

    ## mzXML output:
    out_file <- paste0(test_folder, "/test_write.mzXML")
    suppressWarnings(
        writeMSData(file = out_file, header = hdr, object = pks,
                    outformat = "mzXML")
    )
    in_file <- openMSfile(out_file, backend = "pwiz")
    hdr_2 <- header(in_file)
    pks_2 <- peaks(in_file)
    mzR::close(in_file)
    rt_col <- which(colnames(hdr) == "retentionTime")
    checkEquals(hdr[, rt_col], hdr_2[, rt_col], tolerance = 0.01)
    checkEquals(pks, pks_2)
    hdr_mod <- hdr
    hdr_mod$acquisitionNum <- as.integer(factor(hdr_mod$acquisitionNum))
    hdr_mod$precursorScanNum <- as.integer(factor(hdr_mod$precursorScanNum))
    hdr_2$acquisitionNum <- as.integer(factor(hdr_2$acquisitionNum))
    hdr_2$precursorScanNum <- as.integer(factor(hdr_2$precursorScanNum))
    cn <- colnames(hdr_mod)[!(colnames(hdr_mod) %in% c("retentionTime",
                                                       "injectionTime",
                                                       "filterString",
                                                       "spectrumId"))]
    checkEquals(hdr_mod[, cn], hdr_2[, cn])
    
    ## Subset. These checks ensure that the scan - precursor scan are mapped
    ## correctly!
    idx <- c(`1003` = 1, `1006` = 4, `1019` = 17, `1021` = 19, `1024` = 22,
             `1026` = 24)
    ## 1: no precursor
    ## 2: 1 as precursor, 1003
    ## 3: no precursor
    ## 4: no precursor
    ## 5: 4 as precursor, 1021
    ## 6: 4 as precursor, 1021
    hdr_sub <- hdr[idx, ]
    rownames(hdr_sub) <- NULL
    hdr_sub$seqNum <- 1:nrow(hdr_sub)
    pks_sub <- pks[idx]
    fnew <- paste0(test_folder, "test_copyWrite.mzML")
    writeMSData(file = fnew,
                header = hdr_sub, object = pks_sub, backend = "pwiz",
                outformat = "mzml")
    ## Check content is same
    mzml_new <- openMSfile(fnew, backend = "pwiz")
    pks_new <- peaks(mzml_new)
    hdr_new <- header(mzml_new)
    ii_new <- mzR::instrumentInfo(mzml_new)
    mzR::close(mzml_new)
    checkEquals(pks_new, pks_sub)
    ## I don't quite understand that, but the acquisitionNum and the
    ## precursorScanNum are different while the spectrumId is the same.
    ## Still, check that the precursorScanNum is what we expect:
    checkEquals(hdr_new$precursorScanNum, c(0, 1, 0, 0, 4, 4))
    hdr_new$acquisitionNum <- as.integer(factor(hdr_new$acquisitionNum))
    hdr_new$precursorScanNum <- as.integer(factor(hdr_new$precursorScanNum))
    hdr_sub$acquisitionNum <- as.integer(factor(hdr_sub$acquisitionNum))
    hdr_sub$precursorScanNum <- as.integer(factor(hdr_sub$precursorScanNum))
    checkEquals(hdr_new, hdr_sub)  ## polarity is OK here

    ## Subset with mzXML
    hdr_sub <- hdr[idx, ]
    rownames(hdr_sub) <- NULL
    hdr_sub$seqNum <- 1:nrow(hdr_sub)
    pks_sub <- pks[idx]
    fnew <- paste0(test_folder, "test_copyWrite.mzXML")
    suppressWarnings(
        writeMSData(file = fnew, header = hdr_sub, object = pks_sub,
                    backend = "pwiz", outformat = "mzxml")
    )
    ## Check content is same
    mzml_new <- openMSfile(fnew, backend = "pwiz")
    pks_new <- peaks(mzml_new)
    hdr_new <- header(mzml_new)
    ii_new <- mzR::instrumentInfo(mzml_new)
    mzR::close(mzml_new)
    checkEquals(pks_new, pks_sub)
    checkEquals(hdr_new$precursorScanNum, c(0, 1, 0, 0, 4, 4))
    rownames(hdr_sub) <- NULL
    rownames(hdr_new) <- NULL
    hdr_sub$acquisitionNum <- as.integer(factor(hdr_sub$acquisitionNum))
    hdr_sub$precursorScanNum <- as.integer(factor(hdr_sub$precursorScanNum))
    hdr_new$acquisitionNum <- as.integer(factor(hdr_new$acquisitionNum))
    hdr_new$precursorScanNum <- as.integer(factor(hdr_new$precursorScanNum))
    rt_col <- colnames(hdr_sub) == "retentionTime"
    checkEquals(hdr_sub[, rt_col], hdr_new[, rt_col], tolerance = 0.01)
    cn <- colnames(hdr_sub)[!(colnames(hdr_sub) %in% c("injectionTime",
                                                       "retentionTime",
                                                       "filterString",
                                                       "spectrumId"))]
    checkEquals(hdr_sub[, cn], hdr_new[, cn])
    
    ## Other mzML:
    test_file <- system.file("microtofq", "MM14.mzML", package = "msdata")
    in_file <- openMSfile(test_file, backend = "pwiz")
    hdr <- header(in_file)
    pks <- peaks(in_file)
    mzR::close(in_file)
    
    ## mzML
    out_file <- paste0(test_folder, "/test_write.mzML")
    writeMSData(file = out_file, header = hdr, object = pks)
    in_file <- openMSfile(out_file, backend = "pwiz")
    hdr_2 <- header(in_file)
    pks_2 <- peaks(in_file)
    mzR::close(in_file)
    checkEquals(hdr, hdr_2)
    checkEquals(pks, pks_2)
    ## validate mzML:
    doc <- XML::xmlInternalTreeParse(out_file)
    res <- XML::xmlSchemaValidate(mzML_xsd_idx, doc)
    checkEquals(res$status, 0)

    ## mzXML output:
    out_file <- paste0(test_folder, "test_write.mzXML")
    writeMSData(file = out_file, header = hdr, object = pks,
                outformat = "mzXML")
    in_file <- openMSfile(out_file, backend = "pwiz")
    hdr_2 <- header(in_file)
    pks_2 <- peaks(in_file)
    mzR::close(in_file)
    checkEquals(pks, pks_2)
    checkEquals(hdr[, colnames(hdr_2) != "spectrumId"],
                hdr_2[, colnames(hdr_2) != "spectrumId"])

    ## mzData:
    test_file <- system.file("iontrap", "extracted.mzData", package = "msdata")
    in_file <- openMSfile(test_file, backend = "Ramp")
    hdr <- header(in_file)
    pks <- peaks(in_file)
    mzR::close(in_file)

    ## mzML
    out_file <- paste0(test_folder, "/test_write.mzML")
    writeMSData(file = out_file, header = hdr, object = pks)
    in_file <- openMSfile(out_file, backend = "pwiz")
    hdr_2 <- header(in_file)
    pks_2 <- peaks(in_file)
    mzR::close(in_file)
    checkEquals(hdr, hdr_2)
    checkEquals(pks, pks_2)
    ## validate mzML:
    doc <- XML::xmlInternalTreeParse(out_file)
    res <- XML::xmlSchemaValidate(mzML_xsd_idx, doc)
    checkEquals(res$status, 0)

    ## mzXML output:
    out_file <- paste0(test_folder, "test_write.mzXML")
    writeMSData(file = out_file, header = hdr, object = pks,
                outformat = "mzXML")
    in_file <- openMSfile(out_file, backend = "pwiz")
    hdr_2 <- header(in_file)
    pks_2 <- peaks(in_file)
    mzR::close(in_file)
    checkEquals(pks, pks_2)
    hdr$centroided <- FALSE
    checkEquals(hdr, hdr_2)
}
