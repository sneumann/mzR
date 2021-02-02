setMethod("length",
          signature=c("mzRnetCDF"),
          function(x) {
              return(netCDFVarLen(x@backend, var="scan_number"))
          })

setMethod("peaks", "mzRnetCDF",
          function(object, scans) {
              if (missing(scans))
                  scans <- 1:length(object)
              rawdata <- netCDFRawData(object@backend)

              if (length(scans) == 1) {
                  idx <- seq(rawdata$scanindex[scans] + 1,
                             min(rawdata$scanindex[scans + 1],
                                 length(rawdata$mz), na.rm = TRUE))
                  return(cbind(mz = rawdata$mz[idx],
                               intensity = rawdata$intensity[idx]))
              } else {
                  return(sapply(scans, function(x) {
                      idx <- seq(rawdata$scanindex[x] + 1,
                                 min(rawdata$scanindex[x + 1],
                                     length(rawdata$mz), na.rm = TRUE))
                      cbind(mz = rawdata$mz[idx],
                            intensity = rawdata$intensity[idx])
                  }, simplify = FALSE))
              }
          })

setMethod("spectra", "mzRnetCDF",
          function(object, scans) peaks(object, scans))

## setMethod("peaksCount",
##           signature=c("mzRnetCDF","numeric"),
##           function(object,scans) {
##             if (length(scans)==1) {
##               return(object@backend$getPeakList(scans)$peaksCount)
##             } else {
##               return(sapply(scans,function(x) object@backend$getPeakList(x)$peaksCount))
##             }
##           })

## setMethod("peaksCount",
##           signature=c("mzRnetCDF","missing"),
##           function(object) {
##             n <- length(object)
##             return(peaksCount(object,1:n))
##           })

setMethod("header",
          signature=c("mzRnetCDF","missing"),
          function(object) return(header(object, 1:length(object))))

setMethod("header", c("mzRnetCDF", "numeric"), function(object, scans) {
    ls <- length(scans)
    empty_val <- rep(-1, ls)
    na_real <- rep(NA_real_, ls)
    result <- data.frame(
        seqNum=scans,
        acquisitionNum=scans,
        msLevel=rep(1, length(scans)),
        polarity = empty_val,
        peaksCount=rep(1, length(scans)),
        totIonCurrent=netCDFVarDouble(object@backend, "total_intensity")[scans],
        retentionTime=netCDFVarDouble(object@backend, "scan_acquisition_time")[scans],
        basePeakMZ = empty_val,
        basePeakIntensity = empty_val,
        collisionEnergy = empty_val,
        ionisationEnergy = empty_val,
        lowMZ = empty_val,
        highMZ = empty_val,
        precursorScanNum = empty_val,
        precursorMZ = empty_val,
        precursorCharge = empty_val,
        precursorIntensity = empty_val,
        mergedScan = empty_val,
        mergedResultScanNum = empty_val,
        mergedResultStartScanNum = empty_val,
        mergedResultEndScanNum = empty_val,
        injectionTime = empty_val,
        filterString = NA_character_,
        spectrumId = paste0("scan=", scans),
        centroided = NA,
        ionMobilityDriftTime = empty_val,
        isolationWindowTargetMZ = na_real,
        isolationWindowLowerOffset = na_real,
        isolationWindowUpperOffset = na_real,
        scanWindowLowerLimit = na_real,
        scanWindowUpperLimit = na_real)
    return(result)
})

setMethod("close", 
          signature="mzRnetCDF",
          function(con,...) {
              if (validObject(con))
                  netCDFClose(con@backend)
              con@backend$id <- NULL
              invisible(TRUE)} )

setMethod("isInitialized", 
          signature="mzRnetCDF",
          function(object) return(class(object@backend) == "ncdf4" && validObject(object)))

setMethod("runInfo",
           signature="mzRnetCDF",
          function(object) 
              return(netCDFRunInfo(object@backend)))


setMethod("instrumentInfo",
          signature="mzRnetCDF",
          function(object) 
              return(netCDFInstrumentInfo(object@backend)))


setMethod("manufacturer",
          signature="mzRnetCDF",
          function(object) {
            info <- instrumentInfo(object)           
            return(info$manufacturer)
          })

setMethod("model",
          signature="mzRnetCDF",
          function(object) {
            info <- instrumentInfo(object)           
            return(info$model)
          })

setMethod("ionisation",
          signature="mzRnetCDF",
          function(object) {
            info <- instrumentInfo(object)           
            return(info$ionisation)
          })

setMethod("analyzer",
          signature="mzRnetCDF",
          function(object) {
            info <- instrumentInfo(object)           
            return(info$analyzer)
          })

setMethod("detector",
          signature="mzRnetCDF",
          function(object) {
            info <- instrumentInfo(object)           
            return(info$detector)
          })

setMethod("show",
          signature="mzRnetCDF",
          function(object) {
            if (!isInitialized(object)) {
              cat("Your object's netCDF slot is not initialized.\n")
            } else {
              filename <- fileName(object)
              cat("Mass Spectrometry file handle.\n")
              cat("Filename: ", filename, "\n")
              cat("Number of scans: ", length(object), "\n")
            }
            invisible(NULL)
          })

setMethod("chromatograms", "mzRnetCDF", function(object, chrom)
    chromatogram(object, chrom))
setMethod("chromatogram", "mzRnetCDF", function(object, chrom) {
    warning("The mzRnetCdf backend does not support chromatographic data")
    .empty_chromatogram()
})
setMethod("chromatogramHeader", "mzRnetCDF", function(object, chrom) {
    warning("The mzRnetCdf backend does not support chromatographic data")
    .empty_chromatogram_header()
})
