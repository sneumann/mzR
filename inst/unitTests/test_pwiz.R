test_mzXML <- function() {
    file <- system.file("threonine", "threonine_i2_e35_pH_tree.mzXML", package = "msdata")
    mzxml <- openMSfile(file, backend="pwiz")
    checkTrue(class(mzxml)=="mzRpwiz")

    show(mzxml)
    length(mzxml)
    runInfo(mzxml)
    instrumentInfo(mzxml)
    peaks(mzxml)
    peaks(mzxml,1)
    peaks(mzxml, 2:3)
    peaksCount(mzxml)
    hdr <- header(mzxml)
    checkTrue(any(colnames(hdr) == "spectrumId"))
    header(mzxml,1)
    header(mzxml, 2:3)
    fileName(mzxml)
    close(mzxml) 
}

test_mzML <- function() {
    file <- system.file("microtofq", "MM14.mzML", package = "msdata")
    mzml <- openMSfile(file, backend="pwiz")
    checkTrue(class(mzml)=="mzRpwiz")
    show(mzml)
    length(mzml)
    runInfo(mzml)
    instrumentInfo(mzml)
    pks <- peaks(mzml)
    pks <- peaks(mzml,1)
    pks <- peaks(mzml,2:3)
    peaksCount(mzml)
    hdr <- header(mzml)
    checkTrue(any(colnames(hdr) == "spectrumId"))
    checkEquals(hdr$spectrumId, paste0("spectrum=", hdr$acquisitionNum))
    header(mzml,1)
    header(mzml,2:3)

    checkTrue(ncol(header(mzml))>4)
    checkTrue(length(header(mzml,1))>4)
    checkTrue(ncol(header(mzml,2:3))>4)

    ## Check polarity reporting
    checkTrue(all(header(mzml)$polarity==1))

    fileName(mzml)
    close(mzml)
}

## Test the new implementation of the getScanHeaderInfo
test_getScanHeaderInfo <- function() {
    ## Compare the data returned from the new function with the one returned
    ## by the Ramp backend.
    file <- system.file("microtofq", "MM14.mzML", package = "msdata")
    mzml <- openMSfile(file, backend = "pwiz")
    ramp <- openMSfile(file, backend = "Ramp")
    ## Read single scan header.
    scan_3 <- header(mzml, scans = 3)
    cn <- names(scan_3)
    cn <- cn[cn != "spectrumId"]
    scan_3_ramp <- header(ramp, scans = 3)
    ## Ramp does not read polarity
    scan_3$polarity <- 0
    checkEquals(scan_3[cn], scan_3_ramp[cn])
    
    ## Read all scan header
    all_scans <- header(mzml)
    all_scans_ramp <- header(ramp)
    all_scans$polarity <- 0
    checkEquals(all_scans[, cn], all_scans_ramp[, cn])
    
    ## passing the index of all scan headers should return the same
    all_scans_2 <- header(mzml, scans = 1:nrow(all_scans))
    all_scans_ramp_2 <- header(ramp, scans = 1:nrow(all_scans))
    all_scans_2$polarity <- 0
    checkEquals(all_scans, all_scans_2)
    checkEquals(as.list(all_scans[3, ]), scan_3)
    checkEquals(all_scans_2[, cn], all_scans_ramp_2[, cn])

    ## Some selected scans
    scan_3 <- header(mzml, scans = c(3, 1, 14))
    scan_3_ramp <- header(ramp, scans = c(3, 1, 14))
    ## Ramp does not read polarity
    scan_3$polarity <- 0
    checkEquals(scan_3[, cn], scan_3_ramp[, cn])

    close(mzml)
    close(ramp)

    ## The same for an mzXML file:
    file <- system.file("threonine", "threonine_i2_e35_pH_tree.mzXML",
                        package = "msdata")
    mzml <- openMSfile(file, backend = "pwiz")
    ramp <- openMSfile(file, backend = "Ramp")
    ## Read single scan header.
    scan_3 <- header(mzml, scans = 3)
    scan_3_ramp <- header(ramp, scans = 3)
    checkEquals(scan_3[cn], scan_3_ramp[cn])
    
    ## Read all scan header
    all_scans <- header(mzml)
    all_scans_ramp <- header(ramp)
    ## Ramp unable to read precursorScanNum from an mzXML file.
    all_scans$precursorScanNum <- 0
    checkEquals(all_scans[, cn], all_scans_ramp[, cn])
    
    ## passing the index of all scan headers should return the same
    all_scans_2 <- header(mzml, scans = 1:nrow(all_scans))
    all_scans_ramp_2 <- header(ramp, scans = 1:nrow(all_scans))
    all_scans_2$precursorScanNum <- 0
    checkEquals(all_scans, all_scans_2)
    checkEquals(as.list(all_scans[3, cn]), scan_3[cn])
    checkEquals(all_scans_2[, cn], all_scans_ramp_2[, cn])

    ## Some selected scans
    scan_3 <- header(mzml, scans = c(3, 1, 14))
    scan_3_ramp <- header(ramp, scans = c(3, 1, 14))
    scan_3$precursorScanNum <- 0
    checkEquals(scan_3[, cn], scan_3_ramp[, cn])

    close(mzml)
    close(ramp)

    ## Again for an MSn mzml file.
    file <- msdata::proteomics(full.names = TRUE,
                               pattern = "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz")
    mzml <- openMSfile(file, backend = "pwiz")
    ramp <- openMSfile(file, backend = "Ramp")
    ## Read single scan header.
    scan_3 <- header(mzml, scans = 3)
    scan_3_ramp <- header(ramp, scans = 3)
    ## Ramp does not read polarity, injectionTime or filterString
    scan_3$polarity <- 0
    scan_3$injectionTime <- 0
    scan_3$filterString <- NA_character_
    checkEquals(scan_3[cn], scan_3_ramp[cn])
    
    ## Read all scan header
    all_scans <- header(mzml)
    all_scans_ramp <- header(ramp)
    ## Ramp does not read polarity, injectionTime or filterString
    all_scans$polarity <- 0
    all_scans$injectionTime <- 0
    all_scans$filterString <- NA_character_
    checkEquals(all_scans[, cn], all_scans_ramp[, cn])
    
    ## passing the index of all scan headers should return the same
    all_scans_2 <- header(mzml, scans = 1:nrow(all_scans))
    all_scans_ramp_2 <- header(ramp, scans = 1:nrow(all_scans))
    ## Ramp does not read polarity, injectionTime or filterString
    all_scans_2$polarity <- 0
    all_scans_2$injectionTime <- 0
    all_scans_2$filterString <- NA_character_
    checkEquals(all_scans[, cn], all_scans_2[, cn])
    checkEquals(as.list(all_scans[3, cn]), scan_3[cn])
    checkEquals(all_scans_2[, cn], all_scans_ramp_2[, cn])

    ## Some selected scans
    scan_3 <- header(mzml, scans = c(3, 1, 14))
    scan_3_ramp <- header(ramp, scans = c(3, 1, 14))
    ## Ramp does not read polarity, injectionTime or filterString
    scan_3$polarity <- 0
    scan_3$injectionTime <- 0
    scan_3$filterString <- NA_character_
    checkEquals(scan_3[, cn], scan_3_ramp[, cn])

    close(mzml)
    close(ramp)
}
