#' @title Simple function to check the content of the header argument
#'
#' @description Checks the \code{header} input parameter for the presence of
#'     all required columns before passing it to the C++ routines.
#'
#' @param x a \code{data.frame} in the format as returned by \code{mzR::header}.
#'
#' @author Johannes Rainer
#'
#' @return \code{TRUE} if \code{x} is in the correct format and a
#'     \code{character} with the error message otherwise.
#' 
#' @noRd
.validHeader <- function(x) {
    req_cols <- c(seqNum = "numeric",
                  acquisitionNum = "numeric",
                  msLevel = "numeric",
                  polarity = "numeric",
                  peaksCount = "numeric",
                  totIonCurrent = "numeric",
                  retentionTime = "numeric",
                  basePeakMZ = "numeric",
                  basePeakIntensity = "numeric",
                  collisionEnergy = "numeric",
                  ionisationEnergy = "numeric",
                  lowMZ = "numeric",
                  highMZ = "numeric",
                  precursorScanNum = "numeric",
                  precursorMZ = "numeric",
                  precursorCharge = "numeric",
                  precursorIntensity = "numeric",
                  mergedScan = "numeric",
                  mergedResultScanNum = "numeric",
                  mergedResultStartScanNum = "numeric",
                  mergedResultEndScanNum = "numeric"
                  )
    if (!is.data.frame(x))
        return("'x' is supposed to be a data.frame")
    if (!(all(names(req_cols) %in% colnames(x))))
        return(paste0("'x' is missing one or more required columns: ",
                      paste(names(req_cols), collapse = ", ")))
    cn_x <- colnames(x)
    for (i in 1:ncol(x)) {
        if (!is(x[, i], req_cols[cn_x[i]]))
            return(paste0("column ", cn_x[i], " is expected to contain ",
                          req_cols[cn_x[i]], " values but is of type ",
                          class(x[, i])))
    }
    TRUE
}

#' @title Simple function to check the content of a spectrum list
#'
#' @description Checks whether the passed argument is in the expected format,
#'     which is what is returned by the \code{mzR::peaks} method.
#'
#' @param x a \code{list} such as returned by the \code{mzR::peaks}.
#'
#' @return \code{TRUE} if \code{x} is in the correct format and a
#'     \code{character} with the error message otherwise.
#'
#' @noRd
.validSpectrumList <- function(x) {
    if (!is.list(x))
        return("'x' is supposed to be a list")
    is_ok <- unlist(lapply(x, function(z) {
        if (!is.matrix(z))
            return("list element is not a matrix")
        if (!is.numeric(z))
            return("list should contain only numeric matrices")
        if (ncol(z) != 2)
            return("list should contain matrices with two columns")
        NULL
    }))
    if (length(is_ok))
        return(is_ok[1])
    TRUE
}
