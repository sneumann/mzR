test_mzXML <- function() {
    file <- system.file("threonine", "threonine_i2_e35_pH_tree.mzXML", package = "msdata")
    mzxml <- openMSfile(file, backend="pwiz")
    checkTrue(class(mzxml)=="mzRpwiz")

    show(mzxml)
    length(mzxml)
    runInfo(mzxml)
    instrumentInfo(mzxml)
    checkTrue(is.list(peaks(mzxml)))
    checkTrue(is.matrix(peaks(mzxml,1)))
    pks <- peaks(mzxml, 2:3)
    checkTrue(length(pks) == 2)
    checkEquals(colnames(pks[[1L]]), c("mz", "intensity"))
    peaksCount(mzxml)
    hdr <- header(mzxml)
    checkTrue(any(colnames(hdr) == "spectrumId"))
    checkTrue(all(hdr$centroided))
    checkTrue(any(colnames(hdr) == "scanWindowLowerLimit"))
    checkTrue(any(colnames(hdr) == "scanWindowUpperLimit"))
    checkTrue(all(!is.na(hdr$scanWindowLowerLimit)))
    checkTrue(all(!is.na(hdr$scanWindowUpperLimit)))
    hdr <- header(mzxml,1)
    checkTrue(is.list(hdr))
    hdr <- header(mzxml, 2:3)
    checkTrue(is.data.frame(hdr))
    checkTrue(nrow(hdr) == 2)
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
    checkEquals(colnames(pks[[1L]]), c("mz", "intensity"))
    peaksCount(mzml)
    hdr <- header(mzml)
    checkTrue(any(colnames(hdr) == "spectrumId"))
    checkTrue(all(hdr$centroided))
    checkEquals(hdr$spectrumId, paste0("spectrum=", hdr$acquisitionNum))
    checkTrue(any(colnames(hdr) == "scanWindowLowerLimit"))
    checkTrue(any(colnames(hdr) == "scanWindowUpperLimit"))
    checkTrue(all(!is.na(hdr$scanWindowLowerLimit)))
    checkTrue(all(!is.na(hdr$scanWindowUpperLimit)))
    hdr <- header(mzml,1)
    checkTrue(is.list(hdr))
    hdr <- header(mzml, 2:3)
    checkTrue(is.data.frame(hdr))
    checkTrue(nrow(hdr) == 2)

    checkTrue(ncol(header(mzml))>4)
    checkTrue(length(header(mzml,1))>4)
    checkTrue(ncol(header(mzml,2:3))>4)

    ## Check polarity reporting
    checkTrue(all(header(mzml)$polarity==1))

    fileName(mzml)
    close(mzml)

    file <- system.file("proteomics", "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz",
                        package = "msdata")
    mzml <- openMSfile(file, backend="pwiz")
    checkTrue(class(mzml)=="mzRpwiz")
    hdr <- header(mzml)
    checkTrue(hdr$centroided[2])
    checkTrue(!hdr$centroided[1])
    close(mzml)
}

## Test the new implementation of the getScanHeaderInfo
test_getScanHeaderInfo <- function() {
    ## Compare the data returned from the new function with the one returned
    ## by the Ramp backend.
    file <- system.file("microtofq", "MM14.mzML", package = "msdata")
    mzml <- openMSfile(file, backend = "pwiz")
    ## Read single scan header.
    scan_3 <- header(mzml, scans = 3)
    cn <- names(scan_3)
    cn <- cn[!(cn %in% c("spectrumId", "scanWindowLowerLimit",
                         "scanWindowUpperLimit", "centroided"))]
    scan_3 <- scan_3[cn]
    ## Special case for columns that have an NA reported: they might have a
    ## 0 or NA in Ramp.
    nas <- is.na(scan_3)
    ## Ramp does not read polarity
    scan_3$polarity <- 0

    ## Read all scan header
    all_scans <- header(mzml)
    all_scans <- all_scans[, cn]
    all_scans$polarity <- 0
    nas <- vapply(all_scans, function(z) all(is.na(z)), logical(1))
        
    ## passing the index of all scan headers should return the same
    all_scans <- header(mzml, scans = 1:nrow(all_scans))
    all_scans <- all_scans[, cn]
    all_scans$polarity <- 0
    nas <- vapply(all_scans, function(z) all(is.na(z)), logical(1))
    
    ## Some selected scans
    all_scans <- header(mzml, scans = c(3, 1, 14))
    all_scans <- all_scans[, cn]
    all_scans$polarity <- 0
    nas <- vapply(all_scans, function(z) all(is.na(z)), logical(1))

    close(mzml)

    ## The same for an mzXML file:
    file <- system.file("threonine", "threonine_i2_e35_pH_tree.mzXML",
                        package = "msdata")
    mzml <- openMSfile(file, backend = "pwiz")
    ## Read single scan header.
    scan_3 <- header(mzml, scans = 3)
    scan_3 <- scan_3[cn]
    nas <- vapply(scan_3, function(z) all(is.na(z)), logical(1))
    
    ## Read all scan header
    all_scans <- header(mzml)
    all_scans <- all_scans[, cn]
    
    ## Ramp unable to read precursorScanNum from an mzXML file.
    all_scans$precursorScanNum <- 0
    ## Replace all NA in all_scans with 0
    all_scans <- data.frame(lapply(all_scans, function(z) {
        z[is.na(z)] <- 0
        z
    }))
    
    ## passing the index of all scan headers should return the same
    all_scans <- header(mzml, scans = 1:nrow(all_scans))
    all_scans$precursorScanNum <- 0
    all_scans <- all_scans[, cn]
    ## Replace all NA in all_scans with 0
    all_scans <- data.frame(lapply(all_scans, function(z) {
        z[is.na(z)] <- 0
        z
    }))

    ## Some selected scans
    all_scans <- header(mzml, scans = c(3, 1, 14))
    all_scans$precursorScanNum <- 0
    all_scans <- all_scans[, cn]
    ## Replace all NA in all_scans with 0
    all_scans <- data.frame(lapply(all_scans, function(z) {
        z[is.na(z)] <- 0
        z
    }))

    close(mzml)

    ## Again for an MSn mzml file.
    file <- msdata::proteomics(full.names = TRUE,
                               pattern = "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01.mzML.gz")
    mzml <- openMSfile(file, backend = "pwiz")
    ## Read single scan header.
    scan_3 <- header(mzml, scans = 3)
    ## Ramp does not read polarity, injectionTime or filterString
    cn <- names(scan_3)
    cn <- cn[!(cn %in% c("spectrumId", "scanWindowLowerLimit",
                         "scanWindowUpperLimit", "centroided",
                         "polarity", "filterString",
                         "isolationWindowTargetMZ",
                         "isolationWindowLowerOffset",
                         "isolationWindowUpperOffset",
                         "injectionTime"))]
    scan_3 <- scan_3[cn]
    scan_3[is.na(scan_3)] <- 0 
        
    ## Read all scan header
    all_scans <- header(mzml)[, cn]
    ## Replace all NA in all_scans with 0
    all_scans <- data.frame(lapply(all_scans, function(z) {
        z[is.na(z)] <- 0
        z
    }))
    
    ## passing the index of all scan headers should return the same
    all_scans_2 <- header(mzml, scans = 1:nrow(all_scans))[, cn]
    ## Replace all NA in all_scans with 0
    all_scans_2 <- data.frame(lapply(all_scans_2, function(z) {
        z[is.na(z)] <- 0
        z
    }))
    checkEquals(all_scans, all_scans_2)
    
    ## Some selected scans
    all_scans <- header(mzml, scans = c(3, 1, 14))[, cn]
    ## Replace all NA in all_scans with 0
    all_scans <- data.frame(lapply(all_scans, function(z) {
        z[is.na(z)] <- 0
        z
    }))

    close(mzml)
}
