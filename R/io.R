openMSfile <- function(filename,
                       backend=c("Ramp", "netCDF"),
                       verbose = FALSE) {
    if (!file.exists(filename))
        stop("File ",filename," not found.\n")
    filename <- path.expand(filename)
    if (missing(backend)) {
        ## Guess from file extension
        if (grepl('\\.cdf$', filename,
                  ignore.case = TRUE, perl = TRUE)) {
            backend <- "netCDF"
        } else {
            ## so far everything else is handled by Ramp
            backend <- "Ramp"
        }    
    }
    
    if (tolower(backend) == "ramp") {
        rampModule <- new( Ramp ) 
        rampModule$open(filename, declaredOnly = TRUE)
        if (!rampModule$OK()) {
            stop("Unable to create valid cRamp object.")
        }
        return(new("mzRramp",
                   backend=rampModule,
                   fileName=filename))
    } else if (tolower(backend) == "netcdf") { 
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
    } else {
        stop("No valid backend", backend )
    }  
}



