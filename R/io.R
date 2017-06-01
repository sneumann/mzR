openMSfile <- function(filename,
                       backend=c("Ramp", "pwiz", "netCDF"),
                       verbose = FALSE) {
    if (!file.exists(filename))
        stop("File ",filename," not found.\n")
    filename <- path.expand(filename)
    backend <- match.arg(backend)
    
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
        pwizModule <- new( Pwiz ) 
        pwizModule$open(filename)
        return(new("mzRpwiz",
                   backend=pwizModule,
                   fileName=filename))        
    } else {
        stop("No valid backend", backend)
    }  
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

#' @title Export spectrum data to an MS file
#'
#' @description \code{writeSpectrumList} exports the spectrum data provided with
#'     parameters \code{header} and \code{data} to an MS file.
#'
#' @param filename \code{character(1)} with the name of the file that should be
#'     written.
#'
#' @param header \code{data.frame} with the header data for the spectra. Similar
#'     content as the one returned by the \code{\link{header}} method.
#'
#' @param data \code{list} containing for each spectrum one \code{matrix} with
#'     columns \code{mz} (first column) and \code{intensity} (second column).
#'     See also \code{\link{peaks}} for the method that reads such data from
#'     an MS file.
#'
#' @param backend \code{character(1)} defining the backend that should be used
#'     for writing.
#'
#' @param outformat \code{character(1)} the format of the output file.
#'
#' @param rtime_seconds \code{logical(1)} whether the retention time is provided
#'     in seconds or minutes.
#' 
#' @author Johannes Rainer
writeSpectrumList <- function(filename, backend = "pwiz",
                              outformat = c("mzml", "mgf", "mzxml"), header,
                              data, rtime_seconds = TRUE) {
    backend <- match.arg(backend)
    outformat <- match.arg(tolower(outformat))
    is_ok <- .validHeader(header)
    if (is(is_ok, "character"))
        stop("Error checking parameter 'header': ", is_ok)
    is_ok <- .validSpectrumList(data)
    if (is(is_ok, "character"))
        stop("Error checking parameter 'data'. First error was: ", is_ok)
    ## Ideally, do some checking here.
    if (backend == "pwiz") {
        pwizModule <- new(Pwiz)
        pwizModule$writeSpectrumList(filename, outformat, header, data,
                                     rtime_seconds)
    }
}

#' @title Write spectra data to a MS file
#'
#' @description Copy general information from the originating MS file and
#'     write this, along with the provided spectra data, to a new file. The
#'     expected workflow is the following: data is first loaded from an MS file,
#'     processed in R and then saved again to an MS file.
#'
#' @note This function does not allow to write new MS files with new content.
#'
#' @param filename \code{character(1)} with the name of the file that should be
#'     written.
#'
#' @param originalFile \code{character(1)} with the name of the original file
#'     from which the spectrum data was first read.
#'
#' @param header \code{data.frame} with the header data for the spectra. Similar
#'     content as the one returned by the \code{\link{header}} method.
#'
#' @param data \code{list} containing for each spectrum one \code{matrix} with
#'     columns \code{mz} (first column) and \code{intensity} (second column).
#'     See also \code{\link{peaks}} for the method that reads such data from
#'     an MS file.
#'
#' @param backend \code{character(1)} defining the backend that should be used
#'     for writing.
#'
#' @param outformat \code{character(1)} the format of the output file.
#'
#' @param rtime_seconds \code{logical(1)} whether the retention time is provided
#'     in seconds or minutes.
#' 
#' @author Johannes Rainer
copyWriteMSfile <- function(filename, originalFile, header, data,
                            backend = "pwiz",
                            outformat = "mzml",
                            rtime_seconds = TRUE) {
    backend <- match.arg(backend)
    supp_formats <- c("mzml", "mgf", "mzxml")
    outformat <- match.arg(tolower(outformat), supp_formats)
    if (missing(filename))
        stop("'filename' is a required parameter")
    if (missing(originalFile))
        stop("'originalFile' is a required parameter")
    if (missing(header) | missing(data))
        stop("'header' and 'data' are required")
    if (!file.exists(originalFile))
        stop("Original file ", originalFile, " not found")
    ## Other checks:
    is_ok <- .validHeader(header)
    if (is(is_ok, "character"))
        stop(is_ok)
    is_ok <- .validSpectrumList(data)
    if (is(is_ok, "character"))
        stop("Error checking parameter 'data'. First error was: ", is_ok)
    if (backend == "pwiz") {
        pwizModule <- new(Pwiz)
        pwizModule$copyWriteMSfile(filename, outformat, originalFile,
                                   header, data, rtime_seconds)
    }
}

