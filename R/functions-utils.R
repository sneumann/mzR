#' @title Simple function to validate a 'header' data.frame
#'
#' @description Checks the \code{header} input parameter for the presence of
#'     all required columns before passing it to the C++ routines. The function
#'     in addition ensures that all columns are in the correct order.
#'
#' @param x a \code{data.frame} in the format as returned by \code{mzR::header}.
#'
#' @author Johannes Rainer
#'
#' @return The validated and eventually corrected \code{data.frame} if \code{x}
#'     is in the correct format or a \code{character} with the error message
#'     if not.
#' 
#' @noRd
.validateHeader <- function(x) {
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
                  mergedResultEndScanNum = "numeric",
                  injectionTime = "numeric"
                  )
    if (!is.data.frame(x))
        return("'x' is supposed to be a data.frame")
    if (!(all(names(req_cols) %in% colnames(x))))
        return(paste0("'x' is missing one or more required columns: ",
                      paste(names(req_cols), collapse = ", ")))
    ## Subset and order the columns
    x <- x[, names(req_cols)]
    cn_x <- colnames(x)
    for (i in 1:ncol(x)) {
        if (!is(x[, i], req_cols[cn_x[i]]))
            return(paste0("column ", cn_x[i], " is expected to contain ",
                          req_cols[cn_x[i]], " values but is of type ",
                          class(x[, i])))
    }
    x
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

.check_software_processing <- function(x) {
    if (missing(x))
        return(list())
    if (is.character(x))
        x <- list(x)
    if (is.list(x)) {
        check_element <- function(z) {
            if (!is.character(z))
                stop("Each element in 'software_processing' has to be of type ",
                     "character")
            if (length(z) == 2)
                z <- c(z, "MS:-1")
            if (length(z) < 2)
                stop("Each element in 'software_processing' has to be of ",
                     "length >= 2")
            ## Eventually check that all elements > 2 start with MS.
            z
        }
        x <- lapply(x, check_element)
    } else
        stop("Parameter 'software_processing' has the wrong format")
    x
}
