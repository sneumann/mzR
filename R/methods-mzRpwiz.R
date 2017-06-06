setMethod("get3Dmap", "mzRpwiz",
          function(object, scans, lowMz, highMz, resMz)
          return(object@backend$get3DMap(scans, lowMz, highMz, resMz)))

##setMethod("writeMSfile", "mzRpwiz",
##          function(object, filename, outformat)
##          object@backend$writeMSfile(filename, outformat))

setMethod("length", "mzRpwiz",
          function(x) return(x@backend$getLastScan()))

setMethod("instrumentInfo", "mzRpwiz",
          function(object)
          return(object@backend$getInstrumentInfo()))

setMethod("chromatogramsInfo", "mzRpwiz",
          function(object) {
              .Defunct("chromatogram")
          })


setMethod("manufacturer", "mzRpwiz",
          function(object) {
            info <- instrumentInfo(object)
            return(info$manufacturer)
          })

setMethod("model", "mzRpwiz",
          function(object) {
            info <- instrumentInfo(object)
            return(info$model)
          })

setMethod("ionisation", "mzRpwiz",
          function(object) {
            info <- instrumentInfo(object)
            return(info$ionisation)
          })

setMethod("analyzer", "mzRpwiz",
          function(object) {
            info <- instrumentInfo(object)
            return(info$analyzer)
          })

setMethod("detector", "mzRpwiz",
          function(object) {
            info <- instrumentInfo(object)
            return(info$detector)
          })

setMethod("header", c("mzRpwiz", "missing"),
          function(object) return(object@backend$getAllScanHeaderInfo()))

setMethod("header", c("mzRpwiz", "numeric"),
          function(object, scans) {
              if (length(scans) == 1) {
                  res <- object@backend$getScanHeaderInfo(scans)
                  ## Convert data.frame to list to be conform with old code
                  return(as.list(res))
              } else {
                  return(object@backend$getScanHeaderInfo(scans))
              }
          })

headerFor <- function(object, idx) {
    if (missing(idx))
        stop("Required parameter 'idx' is missing.")
    return(object@backend$getScanHeaderInfoFor(as.integer(idx)))
}

setMethod("peaks", "mzRpwiz",
          function(object, scans) .peaks(object, scans))
setMethod("spectra", "mzRpwiz",
          function(object, scans) .peaks(object, scans))

setMethod("peaksCount", c("mzRpwiz", "numeric"),
          function(object, scans) {
            if (length(scans) == 1) {
              return(object@backend$getPeakList(scans)$peaksCount)
            } else {
                return(sapply(scans,
                              function(x)
                                  object@backend$getPeakList(x)$peaksCount))
            }
          })

setMethod("peaksCount", c("mzRpwiz", "missing"),
          function(object) {
            n <- length(object)
            return(peaksCount(object, 1:n))
          })

setMethod("runInfo", "mzRpwiz",
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

setMethod("softwareInfo", "mzRpwiz",
          function(object) {
            info <- instrumentInfo(object)
            return(info$software)
          })

setMethod("sampleInfo", "mzRpwiz",
          function(object) {
            info <- instrumentInfo(object)
            return(info$sample)
          })

setMethod("sourceInfo", "mzRpwiz",
          function(object) {
            info <- instrumentInfo(object)
            return(info$source)
          })

setMethod("close", "mzRpwiz",
          function(con, ...) {
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

setMethod("isolationWindow", "mzRpwiz",
          function(object, ...) .isolationWindow(fileName(object), ...))

## Chromatograms

nChrom <- function(object) {
          stopifnot(inherits(object, "mzRpwiz"))
          object@backend$getLastChrom()
}

setMethod("tic", "mzRpwiz",
          function(object, ...) {
              if (nChrom(object) < 1)
                  stop("Not chromatogram data available.")
              object@backend$getChromatogramsInfo(0L)
          })

setMethod("chromatograms", "mzRpwiz",
          function(object, chrom) chromatogram(object, chrom))


setMethod("chromatogram", "mzRpwiz",
          function(object, chrom) {
              ## To avoid confusion, the first chromatogram (at index
              ## 0) is indexed at position 1 in R and the last one (at
              ## index nChrom(object) - 1) is indexed at position
              ## nChrom(object).
              n <- nChrom(object)
              all <- FALSE
              if (missing(chrom)) {
                  chrom <- 1:n
                  all <- TRUE
              }
              stopifnot(is.numeric(chrom))
              chrom <- as.integer(chrom)
              if (min(chrom) < 1 | max(chrom) > n)
                  stop("Index out of bound [", 1, ":", n, "].")
              ## Update index to match original indices at the C-level
              chrom <- chrom - 1L
              if (length(chrom) == 1 & !all) {
                  ans <- object@backend$getChromatogramsInfo(chrom)
              } else {
                  ans <- lapply(chrom,
                                function(x)
                                    object@backend$getChromatogramsInfo(x))
              }
              return(ans)
          })
