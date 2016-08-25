setMethod("get3Dmap",
          signature="mzRpwiz",
          function(object,scans,lowMz,highMz,resMz)
          return(object@backend$get3DMap(scans,lowMz,highMz,resMz)))

##setMethod("writeMSfile",
##          signature="mzRpwiz",
##          function(object, filename, outformat)
##          object@backend$writeMSfile(filename, outformat))

setMethod("length",
          signature=c("mzRpwiz"),
          function(x) return(x@backend$getLastScan()))

setMethod("instrumentInfo",
          signature="mzRpwiz",
          function(object)
          return(object@backend$getInstrumentInfo()))

setMethod("chromatogramsInfo",
          signature="mzRpwiz",
          function(object)
          return(object@backend$getChromatogramsInfo()))


setMethod("manufacturer",
          signature="mzRpwiz",
          function(object) {
            info <- instrumentInfo(object)
            return(info$manufacturer)
          })

setMethod("model",
          signature="mzRpwiz",
          function(object) {
            info <- instrumentInfo(object)
            return(info$model)
          })

setMethod("ionisation",
          signature="mzRpwiz",
          function(object) {
            info <- instrumentInfo(object)
            return(info$ionisation)
          })

setMethod("analyzer",
          signature="mzRpwiz",
          function(object) {
            info <- instrumentInfo(object)
            return(info$analyzer)
          })

setMethod("detector",
          signature="mzRpwiz",
          function(object) {
            info <- instrumentInfo(object)
            return(info$detector)
          })

setMethod("header",
          signature=c("mzRpwiz","missing"),
          function(object) return(object@backend$getAllScanHeaderInfo()))

setMethod("header",
          signature=c("mzRpwiz","numeric"),
          function(object, scans) {
            if (length(scans)==1) {
              return(object@backend$getScanHeaderInfo(scans))
            } else {
              return(data.frame(t(sapply(scans,function(x) unlist(object@backend$getScanHeaderInfo(x))))))
            }
          })

setMethod("peaks",
          signature=c("mzRpwiz"),
          function(object, scans) {
              if (missing(scans))
                  scans <- 1:length(object)

              if (length(scans) == 1) {
                  return(object@backend$getPeakList(scans)$peaks)
              } else {
                  return(sapply(scans,
                                function(x) object@backend$getPeakList(x)$peaks,
                                simplify = FALSE))
              }
          })

setMethod("peaksCount",
          signature=c("mzRpwiz","numeric"),
          function(object,scans) {
            if (length(scans)==1) {
              return(object@backend$getPeakList(scans)$peaksCount)
            } else {
              return(sapply(scans,function(x) object@backend$getPeakList(x)$peaksCount))
            }
          })

setMethod("peaksCount",
          signature=c("mzRpwiz","missing"),
          function(object) {
            n <- length(object)
            return(peaksCount(object,1:n))
          })

setMethod("runInfo",
          signature="mzRpwiz",
          function(object) {
            hd <- header(object)
            ll <- list()
            ll$'scanCount' <- length(object)
            ll$'lowMz' <- min(hd$lowMZ)
            ll$'highMz' <- max(hd$highMZ)
            ll$'dStartTime' <- min(hd$retentionTime)
            ll$'dEndTime' <- max(hd$retentionTime)
            ll$'msLevels' <- unique(hd$msLevel)
            return(ll)
          })

setMethod("softwareInfo",
          signature="mzRpwiz",
          function(object) {
            info <- instrumentInfo(object)
            return(info$software)
          })

setMethod("sampleInfo",
          signature="mzRpwiz",
          function(object) {
            info <- instrumentInfo(object)
            return(info$sample)
          })

setMethod("sourceInfo",
          signature="mzRpwiz",
          function(object) {
            info <- instrumentInfo(object)
            return(info$source)
          })

setMethod("close",
          signature = "mzRpwiz",
          function(con,...) {
              con@backend$close()
              invisible(TRUE)
          })




setMethod("show", "mzRpwiz",
          function(object) {
              filename <- fileName(object)
              cat("Mass Spectrometry file handle.\n")
              cat("Filename: ", basename(filename), "\n")
              cat("Number of scans: ", length(object), "\n")
          })

pwiz.version <- function() {
    .Call('mzR_pwiz_version', PACKAGE = 'mzR')
}

