openMSfile <- function(filename,
                       backend = NULL,
                       verbose = FALSE) {
    if (missing(filename))
        stop("'filename' is missing")
    if (length(filename) != 1)
        stop("'filename' has to be of length 1")
    if (!file.exists(filename))
        stop("File ", filename, " not found")
    filename <- path.expand(filename)
    if (is.null(backend))
        backend <- .mzRBackend(filename)
    backend <- match.arg(backend, c("pwiz", "Ramp", "netCDF"))
    
    if (backend == "Ramp") {
        rampModule <- new( Ramp ) 
        rampModule$open(filename, declaredOnly = TRUE)
        if (!rampModule$OK()) {
            stop("Unable to create valid cRamp object.")
        }
        return(new("mzRramp",
                   backend=rampModule,
                   fileName=filename))
    } else if (backend == "netCDF") { 
        if (netCDFIsFile(filename)) {
            ncid <- netCDFOpen(filename)
            if (!is.null(attr(ncid, "errortext"))) {
                stop(attr(ncid, "errortext"))
            }
            return(new("mzRnetCDF",
                       backend=ncid,
                       fileName=filename))
        } else {
            stop("Unable to open netCDF file.")
        }
    } else if (backend == "pwiz") {
        pwizModule <- new(Pwiz)
        tryCatch(pwizModule$open(filename), error = function(e) {
            stop("Can not open file ", filename, "! Original error was: ", e,
                 call. = FALSE)
        })
        return(new("mzRpwiz",
                   backend=pwizModule,
                   fileName=filename))        
    } else {
        stop("No valid backend", backend)
    }  
}

#' @title Define the type of mzR backend to use based on the file name or
#'     content
#'
#' @description Simple helper to define the mzR backend that should/can be used
#'     to read the file.
#'
#' @param x \code{character(1)} representing the file name.
#'
#' @return A \code{character(1)} with the name of the backend (either
#'     \code{"netCDF"}, \code{"Ramp"} or \code{"pwiz"}.
#'
#' @author Johannes Rainer, Sebastian Gibb
#'
#' @noRd
.mzRBackend <- function(x = character()) {
    if (length(x) != 1)
        stop("parameter 'x' has to be of length 1")
    ## Use if/else conditions based on a suggestion from sgibb to avoid loops.
    if (grepl("\\.mzml($|\\.)|\\.mzxml($|\\.)", x, ignore.case = TRUE)) {
        return("pwiz")
    } else if (grepl("\\.mzdata($|\\.)", x, ignore.case = TRUE)) {
        return("Ramp")
    } else if (grepl("\\.cdf($|\\.)|\\.nc($|\\.)", x, ignore.case = TRUE)) {
        return("netCDF")
    } else
        suppressWarnings(.mzRBackendFromContent(x))
}

#' Determine the backend from the (first few lines of the) file content.
#' 
#' @author Johannes Rainer
#'
#' @noRd
.mzRBackendFromContent <- function(x = character()) {
    if (length(x) != 1)
        stop("parameter 'x' has to be of length 1")
    suppressWarnings(
        first_lines <- readLines(x, n = 4)
    )

    if (any(grepl("<mz[X]?ML", first_lines))) {
        return("pwiz")
    } else if (any(grepl("<mzData", first_lines))) {
        return("Ramp")
    } else if (rawToChar(readBin(x, raw(), n = 10)[2:4]) == "HDF") {
        return("pwiz")
    } else if (substr(readBin(x, character(), n = 1), 1, 3) == "CDF") {
        return("netCDF")
    } else
        stop("Could not determine file type for ", x)        
}

openIDfile <- function(filename, verbose = FALSE) {
  if (!file.exists(filename))
    stop("File ",filename," not found.\n")
    
  filename <- path.expand(filename)
  
  identModule <- new(Ident) 
  identModule$open(filename)

    return(new("mzRident",
               backend=identModule,
               fileName=filename))
}

setMethod("writeMSData", signature(object = "list", file = "character"),
          function(object, file, header, backend = "pwiz",
                   outformat = "mzml", rtime_seconds = TRUE,
                   software_processing) {
              backend <- match.arg(backend)
              ## supp_formats <- c("mzml", "mgf", "mzxml")
              supp_formats <- c("mzml", "mzxml")
              outformat <- match.arg(tolower(outformat), supp_formats)
              if (missing(header))
                  stop("'header' is required")
              ## Other checks:
              header <- .validateHeader(header)
              if (is(header, "character"))
                  stop("Error checking parameter 'header': ", header)
              is_ok <- .validSpectrumList(object)
              if (is(is_ok, "character"))
                  stop("Error checking parameter 'object'. First error was: ", is_ok)
              ## Check software_processing:
              software_processing <- .check_software_processing(software_processing)
              ## Add mzR processing:
              mzR <- c("mzR", paste(packageVersion("mzR"), collapse = "."), "MS:-1")
              if (outformat == "mzml")
                  mzR <- c(mzR, "MS:1000544")
              if (outformat == "mzxml")
                  mzR <- c(mzR, "MS:1000545")
              software_processing <- c(software_processing, list(mzR))
              if (backend == "pwiz") {
                  if (outformat == "mzxml" & any(header$injectionTime > 0))
                      warning("mzXML export does not support writing ion injection time")
                  if (outformat == "mzxml" & any(!is.na(header$filterString)))
                      warning("mzXML export does not support writing filter string")
                  pwizModule <- new(Pwiz)
                  pwizModule$writeSpectrumList(file, outformat,
                                               header, object, rtime_seconds,
                                               software_processing)
              }
          })

copyWriteMSData <- function(object, file, original_file, header,
                            backend = "pwiz",
                            outformat = "mzml",
                            rtime_seconds = TRUE,
                            software_processing) {
    backend <- match.arg(backend)
    ## supp_formats <- c("mzml", "mgf", "mzxml")
    supp_formats <- c("mzml", "mzxml")
    outformat <- match.arg(tolower(outformat), supp_formats)
    if (missing(file))
        stop("'file' is a required parameter")
    if (missing(original_file))
        stop("'original_file' is a required parameter")
    if (missing(header) | missing(object))
        stop("'header' and 'object' are required")
    if (!file.exists(original_file))
        stop("Original file ", original_file, " not found")
    ## Other checks:
    header <- .validateHeader(header)
    if (is(header, "character"))
        stop("Error checking parameter 'header': ", header)
    is_ok <- .validSpectrumList(object)
    if (is(is_ok, "character"))
        stop("Error checking parameter 'object'. First error was: ", is_ok)
    ## Check software_processing:
    software_processing <- .check_software_processing(software_processing)
    ## Add mzR processing:
    mzR <- c("mzR", paste(packageVersion("mzR"), collapse = "."), "MS:-1")
    if (outformat == "mzml")
        mzR <- c(mzR, "MS:1000544")
    if (outformat == "mzxml")
        mzR <- c(mzR, "MS:1000545")
    software_processing <- c(software_processing, list(mzR))
    if (backend == "pwiz") {
        if (outformat == "mzxml" & any(header$injectionTime > 0)) {
            warning("mzXML export does not support writing ion injection time")
            header$injectionTime = 0
        }
        if (outformat == "mzxml" & any(!is.na(header$filterString))) {
            warning("mzXML export does not support writing filter string")
            header$filterString <- NA_character_
        }
        pwizModule <- new(Pwiz)
        pwizModule$copyWriteMSfile(file, outformat, original_file,
                                   header, object, rtime_seconds,
                                   software_processing)
    }
}

