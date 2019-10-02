.peaks <- function(object, scans) {
    if (missing(scans))
        scans <- 1:length(object)
    res <- object@backend$getPeakList(scans)
    if (length(res) == 1)
        res[[1]]
    else res
}

.peaks_ramp <- function(object, scans) {
    if (missing(scans))
        scans <- 1:length(object)
    if (length(scans) == 1) {
        object@backend$getPeakList(scans)$peaks
    } else {
        sapply(scans, function(x) object@backend$getPeakList(x)$peaks,
               simplify = FALSE)
    }
}

setMethod("isolationWindow", "character",
          function(object, ...) .isolationWindow(object, ...))

.isolationWindow <- function(x, unique. = TRUE, simplify = TRUE) {
    stopifnot(all(file.exists(x)))
    if (!requireNamespace("XML"))
        stop("Please install the XML package to use this functionality.")
    res <- lapply(x, function(xx) {
        xml <- XML::xmlParse(xx)
        ns <- c(x = "http://psi.hupo.org/ms/mzml")
        path <- c(low = "//x:isolationWindow/x:cvParam[@accession='MS:1000828']/@value",
                  high = "//x:isolationWindow/x:cvParam[@accession='MS:1000829']/@value")
        low <- as.numeric(XML::xpathSApply(xml, path["low"], namespaces = ns))
        high <- as.numeric(XML::xpathSApply(xml, path["high"], namespaces = ns))
        cbind(low, high)
    })
    if (.multipleIsolationWindows(res))
        message("Found multiple isolation windows in an acquisition.")
    if (unique.)
        res <- lapply(res, base::unique)
    if (simplify & length(x) == 1) res <- res[[1]]
    return(res)
}

.multipleIsolationWindows <- function(x) {
    x <- lapply(x, base::unique)
    any(sapply(x, function(xx) nrow(xx) > 1))
}


.hasSpectra <- function(x) {
    close_after <- FALSE
    if (is.character(x) && file.exists(x)) {
        x <- mzR::openMSfile(x)
        ## Ensure we are closing the file later
        close_after <- TRUE
    }
    stopifnot(inherits(x, "mzR"))
    len <- length(x)
    if (close_after)
        close(x)
    return(as.logical(len))
}

#' Create return data for MS backends not supporting chromatographic data. This
#' function is supposed to be called by the chromatogram(s) methods for these
#' backends
#'
#' @author Johannes Rainer
#'
#' @noRd
.empty_chromatogram <- function() {
    list()
}

#' Create return data for MS backends not supporting chromatographic data.
#'
#' @author Johannes Rainer
#'
#' @noRd
.empty_chromatogram_header <- function() {
    cn <- c("chromatogramId", "chromatogramIndex", "polarity",
            "precursorIsolationWindowTargetMZ",
            "precursorIsolationWindowLowerOffset",
            "precursorIsolationWindowUpperOffset",
            "precursorCollisionEnergy",
            "productIsolationWindowTargetMZ",
            "productIsolationWindowLowerOffset",
            "productIsolationWindowUpperOffset")
    data.frame(matrix(nrow = 0, ncol = length(cn),
                      dimnames = list(character(), cn)))
}

.hasChromatograms <- function(x) {
    close_after <- FALSE
    if (is.character(x) && file.exists(x)) {
        x <- mzR::openMSfile(x)
        close_after <- TRUE
    }
    stopifnot(inherits(x, "mzR"))
    hdr <- chromatogramHeader(x)
    if (close_after)
        close(x)
    as.logical(nrow(hdr))
}
