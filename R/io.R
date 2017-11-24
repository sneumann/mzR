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

