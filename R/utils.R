.peaks <- function(object, scans) {
    if (missing(scans))
        scans <- 1:length(object)
    if (length(scans) == 1) {
        return(object@backend$getPeakList(scans)$peaks)
    } else {
        return(sapply(scans,
                      function(x) object@backend$getPeakList(x)$peaks,
                      simplify = FALSE))
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
