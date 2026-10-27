#' Read mass spectrometry raw data
#'
#' Open netCDF, mzXML or mzML files using a format-specific backend. Access
#' spectra, chromatograms, and run and instrument metadata through mzR objects.
#'
#' @edam_operation Data parsing (operation_1812)
#' @edam_input_data Mass spectrometry data (data_2536)
#' @edam_input_format mzML (format_3244)
#' @edam_input_format mzXML (format_3654)
#' @edam_input_format NetCDF (format_3650)
#' @edam_output_data Mass spectrometry data (data_2536)
#' @noRd
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
    backend <- match.arg(backend, c("pwiz", "netCDF"))
    
    if (backend == "netCDF") { 
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
#'     \code{"netCDF"} or \code{"pwiz"}.
#'
#' @author Johannes Rainer, Sebastian Gibb
#'
#' @noRd
.mzRBackend <- function(x = character()) {
    if (length(x) != 1)
        stop("parameter 'x' has to be of length 1")
    ## Use if/else conditions based on a suggestion from sgibb to avoid loops.
    if (grepl("\\.mzml($|\\.)|\\.mzxml($|\\.)", x, ignore.case = TRUE))
        return("pwiz")
    else if (grepl("\\.cdf($|\\.)|\\.nc($|\\.)", x, ignore.case = TRUE))
        return("netCDF")
    else
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
    } else if (rawToChar(readBin(x, raw(), n = 10)[2:4]) == "HDF") {
        return("pwiz")
    } else if (substr(readBin(x, character(), n = 1), 1, 3) == "CDF") {
        return("netCDF")
    } else
        stop("Could not determine file type for ", x)        
}

#' Read peptide identification data
#'
#' Open mzIdentML 1.1 files using the ProteoWizard backend to access
#' peptide-spectrum matches, scores and identification metadata.
#'
#' @edam_operation Data parsing (operation_1812)
#' @edam_input_data Mass spectrometry data (data_2536)
#' @edam_input_format mzIdentML (format_3247)
#' @edam_output_data Mass spectrometry data (data_2536)
#' @noRd
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

#' Write mass spectrometry spectra
#'
#' Export a list of spectrum matrices and a header data frame to an mzML or
#' mzXML file, including software processing information.
#'
#' @edam_operation Data handling (operation_2409)
#' @edam_input_data Mass spectrometry data (data_2536)
#' @edam_output_data Mass spectrometry data (data_2536)
#' @edam_output_format mzML (format_3244)
#' @edam_output_format mzXML (format_3654)
#' @noRd
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
              mzR <- c("mzR", paste(packageVersion("mzR"), collapse = "."),
                       "MS:1002869")
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

#' Copy metadata and write modified mass spectrometry data
#'
#' Write spectra in mzML or mzXML format, preserving general metadata from an
#' original mass spectrometry file.
#'
#' @edam_operation Conversion (operation_3434)
#' @edam_input_data Mass spectrometry data (data_2536)
#' @edam_input_format mzML (format_3244)
#' @edam_input_format mzXML (format_3654)
#' @edam_output_data Mass spectrometry data (data_2536)
#' @edam_output_format mzML (format_3244)
#' @edam_output_format mzXML (format_3654)
#' @noRd
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
    mzR <- c("mzR", paste(packageVersion("mzR"), collapse = "."), "MS:1002869")
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
