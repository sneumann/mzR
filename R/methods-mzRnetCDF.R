## setMethod("get3Dmap",
##           signature="mzRnetCDF",
##           function(object,scans,lowMz,highMz,resMz) 
##           return(object@backend$get3DMap(scans,lowMz,highMz,resMz)))

## setMethod("initializeRamp",
##           signature="mzRnetCDF",
##           function(object) {
##             if (!file.exists(fileName(object)))
##               stop("File ",fileName(object)," not found.\n")
##             object@backend$open(fileName(object), declaredOnly = TRUE)
##             if (isInitialized(object)) invisible(TRUE)
##             else stop("Could not initialize ramp slot.")
##           })

setMethod("length",
          signature=c("mzRnetCDF"),
          function(x) {
            scanindex <- netCDFVarInt(x@backend, "scan_index")
            if (!is.null(attr(scanindex, "errortext")))
              stop("Couldn't read scan indicies from ", x@backend)
            return(length(scanindex))
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

setMethod("scans", "mzRnetCDF",
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

setMethod("header",
          signature=c("mzRnetCDF","numeric"),
          function(object, scans) {
            
            result <- data.frame(seqNum=scans,
                            acquisitionNum=scans,
                            msLevel=rep(1, length(scans)),
                            peaksCount=rep(1, length(scans)),
                            totIonCurrent=netCDFVarDouble(object@backend, "total_intensity")[scans],
                            retentionTime=netCDFVarDouble(object@backend, "scan_acquisition_time")[scans], 
                            basePeakMZ=rep(-1, length(scans)),
                            basePeakIntensity=rep(-1, length(scans)),
                            collisionEnergy=rep(-1, length(scans)),
                            ionisationEnergy=rep(-1, length(scans)),
                            highMZ=rep(-1, length(scans)),
                            precursorScanNum=rep(-1, length(scans)),
                            precursorMZ=rep(-1, length(scans)),
                            precursorCharge=rep(-1, length(scans)),
                            precursorIntensity=rep(-1, length(scans)),
                            mergedScan=rep(-1, length(scans)),
                            mergedResultScanNum=rep(-1, length(scans)),
                            mergedResultStartScanNum=rep(-1, length(scans)),
                            mergedResultEndScanNum=rep(-1, length(scans)))

            return(result)
          })

setMethod("close", 
          signature="mzRnetCDF",
          function(con,...) return( netCDFClose(con@backend) ))

setMethod("isInitialized", 
          signature="mzRnetCDF",
          function(object) return(object@backend > 0))

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

