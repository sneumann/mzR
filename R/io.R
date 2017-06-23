openMSfile <- function(filename,
                       backend=c("pwiz", "Ramp", "netCDF"),
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

#' @title Write MS spectrum data to an MS file
#'
#' @description \code{writeMSData} exports the MS spectrum data provided with
#'     parameters \code{header} and \code{data} to an MS file in mzML or
#'     mzXML format.
#'
#' @param filename \code{character(1)} defining the name of the file.
#'
#' @param header \code{data.frame} with the header data for the spectra. Has to
#'     be in the format as the \code{data.frame} returned by the
#'     \code{\link{header}} method.
#'
#' @param data \code{list} containing for each spectrum one \code{matrix} with
#'     columns \code{mz} (first column) and \code{intensity} (second column).
#'     See also \code{\link{peaks}} for the method that reads such data from
#'     an MS file.
#'
#' @param backend \code{character(1)} defining the backend that should be used
#'     for writing. Currently only \code{"pwiz"} backend is supported.
#'
#' @param outformat \code{character(1)} the format of the output file. One of
#'     \code{"mzml"}, \code{"mzxml"} and \code{"mgf"}.
#'
#' @param rtime_seconds \code{logical(1)} whether the retention time is provided
#'     in seconds or minutes (defaults to \code{TRUE}).
#' 
#' @param software_processing \code{list} of \code{character} vectors (or single
#'     \code{character} vector). Each \code{character} vector providing
#'     information about the software that was used to process the data with
#'     optional additional description of processing steps. The length of each
#'     \code{character} vector has to be >= 3: the first element being the name
#'     of the software, the second string its version and the third element the
#'     MS CV ID of the software (or \code{"MS:-1"} if not known). All additional
#'     elements are optional and represent the MS CV ID of each processing step
#'     performed with the software.
#' 
#' @author Johannes Rainer
writeMSData <- function(filename, header, data, backend = "pwiz",
                        outformat = c("mzml"),
                        rtime_seconds = TRUE,
                        software_processing) {
    backend <- match.arg(backend)
    ## supp_formats <- c("mzml", "mgf", "mzxml")
    supp_formats <- "mzml"
    outformat <- match.arg(tolower(outformat), supp_formats)
    if (missing(filename))
        stop("'filename' is a required parameter")
    if (missing(header) | missing(data))
        stop("'header' and 'data' are required")
    ## Other checks:
    header <- .validateHeader(header)
    if (is(header, "character"))
        stop("Error checking parameter 'header': ", header)
    is_ok <- .validSpectrumList(data)
    if (is(is_ok, "character"))
        stop("Error checking parameter 'data'. First error was: ", is_ok)
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
        pwizModule <- new(Pwiz)
        pwizModule$writeSpectrumList(filename, outformat,
                                     header, data, rtime_seconds,
                                     software_processing)
    }
}

#' @title Write MS spectrum data to a MS file copying metadata from the
#'     originating file
#'
#' @description Copy general information from the originating MS file and
#'     write this, along with the provided spectra data, to a new file. The
#'     expected workflow is the following: data is first loaded from an MS file,
#'     e.g. using \code{\link{peaks}} and \code{\link{header}} methods,
#'     processed in R and then saved again to an MS file providing the
#'     (eventually) manipulated spectra and header data with arguments
#'     \code{header} and \code{data}.
#'
#' @note This function does not allow to write new MS files with new content.
#'     Use the \code{\link{writeMSData}} function for that.
#'
#' @inheritParams writeMSData
#' 
#' @param originalFile \code{character(1)} with the name of the original file
#'     from which the spectrum data was first read.
#'
#' @seealso \code{\link{writeMSData}} for a function to save MS data to a new
#'     mzML or mzXML file.
#' 
#' @author Johannes Rainer
copyWriteMSData <- function(filename, originalFile, header, data,
                            backend = "pwiz",
                            outformat = "mzml",
                            rtime_seconds = TRUE,
                            software_processing) {
    backend <- match.arg(backend)
    ## supp_formats <- c("mzml", "mgf", "mzxml")
    supp_formats <- "mzml"
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
    header <- .validateHeader(header)
    if (is(header, "character"))
        stop("Error checking parameter 'header': ", header)
    is_ok <- .validSpectrumList(data)
    if (is(is_ok, "character"))
        stop("Error checking parameter 'data'. First error was: ", is_ok)
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
        pwizModule <- new(Pwiz)
        pwizModule$copyWriteMSfile(filename, outformat, originalFile,
                                   header, data, rtime_seconds,
                                   software_processing)
    }
}

