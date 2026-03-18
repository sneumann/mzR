test_header_all <- function() {
    ## Check that the header call returns the same columns irrespectively of the
    ## backend. Issue #238
    file <- MsDataHub::ko15.CDF()
    cdf <- openMSfile(file, backend = "netCDF")
    header_cdf <- header(cdf)
    close(cdf)

    file <- MsDataHub::PestMix1_DDA.mzML()
    mzxml <- openMSfile(file, backend="pwiz")
    header_pwiz <- header(mzxml)
    close(mzxml)

    neededCdfHeaders <- c("seqNum", "acquisitionNum", "msLevel", "polarity",
                          "peaksCount", "totIonCurrent", "retentionTime",
                          "basePeakMZ", "basePeakIntensity", "collisionEnergy",
                          "ionisationEnergy", "lowMZ", "highMZ",
                          "precursorScanNum", "precursorMZ", "precursorCharge",
                          "precursorIntensity", "mergedScan",
                          "mergedResultScanNum", "mergedResultStartScanNum",
                          "mergedResultEndScanNum", "injectionTime",
                          "filterString", "spectrumId", "centroided",  
                          "ionMobilityDriftTime", "isolationWindowTargetMZ",
                          "isolationWindowLowerOffset",
                          "isolationWindowUpperOffset", "scanWindowLowerLimit",
                          "scanWindowUpperLimit", "electronBeamEnergy")
    
    checkTrue( all(neededCdfHeaders %in% colnames(header_cdf)) )
    checkTrue( all(neededCdfHeaders %in% colnames(header_pwiz)) )
}
