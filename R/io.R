openMSfile <- function(filename,
                       backend=c("Ramp", "netCDF", "pwiz"),
                       verbose = FALSE) {
  if (!file.exists(filename))
    stop("File ",filename," not found.\n")
  filename <- path.expand(filename)
  
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
##      on.exit(netCDFClose(ncid)) ## Hm, that worked in xcms ?!
      return(new("mzRnetCDF",
                 backend=ncid,
                 fileName=filename))
    } else {
      stop("Unable to open netCDF file.")
    }
  }else if(tolower(backend) == "pwiz"){
      pwizModule <- new( Pwiz ) 
	  pwizModule$open(filename)

    return(new("mzRpwiz",
               backend=pwizModule,
               fileName=filename))
  
  } else {
    stop("No valid backend", backend )
  }  
}



